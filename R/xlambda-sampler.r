#' Joint x-lambda sampler for Poisson and negative binomial linear inverse problems
#'
#' Consider the linear inverse problem y=Ax, where x follows a Poisson or negative binomial distribution with mean lambda = E[x]. This function jointly samples lambda and x. 
#' @param y Matrix of sequence of observed count data vectors; each column is an observation.
#' @param A Model configuration matrix, assumed to be binary.
#  @param lambda.updater Function with required arguments x, lambda (and NB.alpha for negative binomial models) which updates lambda. 
#' @param lambda.ini Initial mean vector for x.
#' @param U Optional matrix the columns of which should be a Markov (sub)-basis.
#' @param Method "MH" for Metropolis-Hastings sampler, "Gibbs" for Gibbs sampler.
#' @param Reorder Should the columns of A be reordered? Defaults to TRUE.
#' @param tune.par Tuning parameter (alpha) controlling variation in fitness values for lattice bases. Defaults to 0.5.
#' @param combine Should extra moves be included combining lattice basis vectors? Defaults to FALSE, but should usually be set to TRUE if A is not unimodular.
#' @param x.order If Reorder=FALSE, x.order can be used to reorder columns of A to match ordering of entries of x. Defaults to NULL when no such reordering is performed.
#' @param x.ini Matrix of initial values for x, with column orderings matching that for y. Default is NULL, when initial values derived through integer programming.
#' @param Model "Poisson" or "NegBin".
#' @param Proposal "NonUnif" or "Unif" (default).
#' @param NB.alpha.ini Initial value for dispersion parameter for negaqtive-binomial distribution. Defaults to 1.
#' @param lambda.additional Optional object to transfer additional information to lambda.updater. 
#' @param other.pars Optional object to provide initial values for other model parameters.
#' @param ndraws Number of iterations to run sampler after burn-in. One iteration comprises cycling through the full basis (possibly augmented by a combined move). Defaults to 10^4.
#' @param burnin Number of iteractions for burn in period. Defaults to 2000, which is usually more than adequate.
#' @param verbose Controls level of detail in recording lattice bases used.
#' @param THIN Thinning parameter for output. Defaults to 1 (no thinning).
#' @return A list with components X (an array, for which X[i,j,k] is the k-th sampled value of the i-th component of the j-th observation of x), LAMBDA (a matrix, each row corresponding to samples for an entry of lambda), NB.ALPHA (a vector of sampler values of NB.alpha, NA if model is Poisson), OTHER.PARS (a matrix, each row corresponding to an additional parameter; zero rows if there are none) and x.order (a vector describing dynamic selection of lattice bases, if verbose=1).
#' @export
#' @examples 
#' data(LondonRoad)
#' lu <- function(x,lambda,NB.alpha=NA,lambda.tuning=1,lambda.additional=NA) { list(lambda=rgamma(length(lambda),shape=x+0.5*LondonRoad$lambda,rate=1.5),other.pars=numeric(0),NB.alpha=NA) }
#' Xlambdasampler(y=LondonRoad$y,A=LondonRoad$A,lambda.updater=lu,lambda.ini=LondonRoad$lambda,Model="Poisson",Method="Gibbs",tune.par=0.5,combine=FALSE)

Xlambdasampler <- function (y, A, lambda.updater, lambda.ini, U=NULL, Method="MH", Reorder=TRUE, tune.par=0.5, combine=FALSE, x.order=NULL, x.ini=NULL, Model="Poisson", Proposal="Unif", NB.alpha.ini=1, lambda.additional=NA, other.pars=numeric(0), ndraws = 10000, burnin = 2000, verbose = 0, THIN = 1) {
	require(lpSolve)
	require(numbers)
	require(extraDistr)
  	if(Model=="NegBin" & NB.alpha.ini<=0) NB.alpha.ini=1
	if(Model=="Uniform") Method <- "Gibbs"

	zero.cols.ind <- which(colSums(A)==0)
	non.zero.cols.ind <- which(colSums(A) > 0)
	zero.cols <- (length(zero.cols.ind) > 0)

	Y <- as.matrix(y)
	r <- ncol(A)
	n <- nrow(A)
	ntime <- ncol(Y)

	if (zero.cols) max.ratio <- max(c(5*lambda.ini,c(Y)))

	if(!is.null(x.ini)) x.ini <- as.matrix(x.ini)
	tol <- 10^{-10}

	if (is.matrix(U)) Reorder=FALSE

	lambda <- lambda.ini
	

	if (Reorder){
		lambda.star <- lambda - (colSums(A)==0)*1e15
		lam.order <- order(lambda.star,decreasing=TRUE)  
		A <- A[,lam.order]
		x.order <- lam.order
		for (i in 2:n){
			while (qr(A[,1:i])$rank < i){
				A <- A[,c(1:(i-1),(i+1):r,i)]
				x.order <- x.order[c(1:(i-1),(i+1):r,i)]
			}
		}
	}


	if (!Reorder) { 
		if(is.null(x.order)) x.order <- 1:r
		A <- A[,x.order]
	}

	lambda <- lambda[x.order]
	
	if (is.matrix(U)){
	 	dA1 <- 1
		tune.par <- -1
	}
	if (!is.matrix(U)){
		A1 <- A[,1:n]
		dA1 <- det(A1)
		A2 <- as.matrix(A[,-c(1:n)])
		A1.inv <- solve(A1)
		C <- A1.inv%*%A2
		U <- rbind(-C,diag(r-n))
	}

	m <- ncol(U) + 1*combine

	X <- array(0, dim=c(r,ntime,ndraws + burnin))
	LAMBDA <- matrix(0, r, ndraws + burnin)
	NB.ALPHA <- numeric(ndraws + burnin)
	OTHER.PARS <- matrix(0,length(other.pars),ndraws + burnin)
	OTHER.PARS[,1] <- other.pars
	if (verbose==1) X.ORDER <- matrix(0, r, ndraws + burnin)

	NB.alpha <- NB.alpha.ini

	if (!is.null(x.ini)) x.ini <- x.ini[x.order,]
	if (is.null(x.ini)){
       		x.ini <- matrix(0, nrow = r, ncol = ntime)
	  	A.nozero <- A[,colSums(A)>0.01]
	  	r.nozero <- ncol(A.nozero)
        	for (tt in 1:ntime) {
            		x.ini[colSums(A)>0.01, tt] <- lp("max", objective.in = rep(1, r.nozero), const.mat = A.nozero, const.dir = rep("=", nrow(A)), const.rhs = c(Y[, tt]), all.int = T)$solution
			x.ini[colSums(A)<0.01 ,tt] <- mean(x.ini[colSums(A)>0.01 ,tt])
        	}
	}
	xx <- x.ini
        X[x.order,,1] <- xx
	LAMBDA[x.order,1] <- lambda
	if (verbose==1) X.ORDER[,1] <- x.order

	for (iter in seq(2, ndraws + burnin)) {
		if (tune.par > 0){
                  	lambda.star <- rnorm(r,mean=lambda,sd=tune.par*lambda)
			ii <- sample(1:n,1)
			swap.indx <- abs(C[ii,])>tol
			if (any(swap.indx)){
				if (sum(swap.indx) > 1) {
					jj <- sample((1:(r-n))[swap.indx],1)
				} else {
					jj <- (1:(r-n))[swap.indx]
				}
				if (lambda.star[ii] <= lambda.star[jj+n]/abs(C[ii,jj])){
				ei <- rep(0,n)
 				ej <- rep(0,r-n)
				ei[ii] <- 1
				ej[jj] <- 1
				dA1 <- round(C[ii,jj]*dA1)
				C <- C - outer(C[,jj]-ei,C[ii,]+ej)/C[ii,jj]
				U <- rbind(-C,diag(r-n))
				x.order[c(ii,jj+n)] <- x.order[c(jj+n,ii)]
				xx[c(ii,jj+n),] <- xx[c(jj+n,ii),]
 				lambda[c(ii,jj+n)] <- lambda[c(jj+n,ii)]
 				}
			}
		}
		for (tt in 1:ntime){
			x <- xx[,tt]
			for (j in 1:m){
				if (j <= ncol(U)) z <- U[,j]
				if (j > ncol(U)){
					delta <- 0
					while(sum(delta)==0) delta <- rpois(ncol(U),lambda=0.5)
					delta <- delta*sample(c(-1,1),size=ncol(U),replace=T)
					z <- U%*%delta
				}
				if (abs(dA1)!=1) { if(is_wholenumber(z)==F) z <- round(z*abs(dA1))/mGCD(round(abs(z*dA1))) }
				ratio <- (x/abs(z))[which(z<0)]
				if (length(ratio)==0) ratio <- max.ratio
				max.move <- floor(min(ratio))
				min.move <- -floor(min((x/abs(z))[which(z>0)]))
				x.min <- x+min.move*z
				x.max <- x+max.move*z
				indx <- 1:(max.move-min.move+1)
				update.indx <- which(z!=0)
				if (max(indx>1)){
					if (Method=="Gibbs"){
						if (Model!="Uniform") x.matrix <- round(t(mapply(seq,from=x.min[update.indx],by=z[update.indx],length.out=max.move-min.move+1)))
						if (Model=="Poisson") {
							log.probs <- colSums(dpois(x.matrix,lambda[update.indx],log=T))
							probs <- exp(log.probs-max(log.probs))
							x <- x.min + sample(indx-1,size=1,prob=probs)*z
						}
						if (Model=="NegBin"){
							log.probs <- colSums(dnbinom(x.matrix,mu=lambda[update.indx],size=lambda[update.indx]/NB.alpha,log=T))
							probs <- exp(log.probs-max(log.probs)) 
							x <- x.min + sample(indx-1,size=1,prob=probs)*z
						}
						if (Model=="Uniform") x <- x.min+sample(indx-1,size=1)*z
					}
					if (Method=="MH"){
						if (Proposal=="Unif") move.length <- sample(min.move:max.move,1)
						if (Proposal=="NonUnif"){
	            					if (Model=="Poisson") {
								aa <- x[n+j]+z[n+j]*min.move
								bb <- x[n+j]+z[n+j]*max.move
								move.length <- (rtpois(1,lambda=lambda[n+j],a=aa-0.5,b=bb)-x[n+j])/z[n+j]
							}
	            					if (Model=="NegBin") move.length <- sample(min.move:max.move,1,prob=dnbinom((x[n+j]+z[n+j]*min.move):(x[n+j]+z[n+j]*max.move),mu=lambda[n+j],size=lambda[n+j]/NB.alpha))  			
							if (Model=="Normal") move.length <- rtnorm(1,lambda[n+j],sd=sqrt((1+NB.alpha)*lambda[n+j]),lower=x[n+j]+z[n+j]*min.move,upper=x[n+j]+z[n+j]*max.move)
						} 
						x.cand <- x + z*move.length
						if (Model=="Poisson"){
							L <- sum(dpois(x[update.indx],lambda[update.indx],log=T))
							L.cand <- sum(dpois(x.cand[update.indx],lambda[update.indx],log=T))
						}
						if (Model=="NegBin"){
							L <- sum(dnbinom(x[update.indx],mu=lambda[update.indx],size=lambda[update.indx]/NB.alpha,log=T))
							L.cand <- sum(dnbinom(x.cand[update.indx],mu=lambda[update.indx],size=lambda[update.indx]/NB.alpha,log=T))
						}
						if (Model=="Normal"){
							L <- sum(dnorm(x[update.indx],mean=lambda[update.indx],sd=sqrt((1+NB.alpha)*lambda[update.indx]),log=T))
							L.cand <- sum(dnorm(x.cand[update.indx],mean=lambda[update.indx],sd=sqrt((1+NB.alpha)*lambda[update.indx]),log=T))
						}
						if (Proposal=="Unif") acc.prob <- exp(L.cand - L)
						if (Proposal=="NonUnif"){
	            					if (Model=="Poisson"){
								q.can <- dpois(x.cand[n+j],lambda[n+j],log=T)  
								q.cur <- dpois(x[n+j],lambda[n+j],log=T)
							}
	            					if (Model=="NegBin"){
								q.can <- dnbinom(x.cand[n+j],mu=lambda[n+j],size=lambda[n+j]/NB.alpha,log=T)  
								q.cur <- dnbinom(x[n+j],mu=lambda[n+j],size=lambda[n+j]/NB.alpha,log=T)
							}
 							if (Model=="Normal") {
								q.can <- dnorm(x.cand[n+j],lambda[n+j],sqrt((1+NB.alpha)*lambda[n+j]),log=T)  
								q.cur <- dnorm(x[n+j],lambda[n+j],sqrt((1+NB.alpha)*lambda[n+j]),log=T)
							} 
			  				acc.prob <- exp(L.cand - L + q.cur - q.can)
						}
						if (is.na(acc.prob)) acc.prob <- 0
						if (runif(1) < acc.prob) x <- x.cand
					}
				}

			}
		X[x.order,tt,iter] <- x
		xx[,tt] <- x
		}
		if (verbose==1) X.ORDER[,iter] <- x.order
		if (Model=="NegBin"){
			updates <- lambda.updater(x=xx[order(x.order),],lambda=lambda[order(x.order)],NB.alpha=NB.alpha,lambda.additional=lambda.additional,other.pars=other.pars)
			lambda <- updates$lambda[x.order]
			NB.alpha <- updates$NB.alpha
			lambda.additional <- updates$lambda.additional
			other.pars <- updates$other.pars
			LAMBDA[,iter] <- updates$lambda
			NB.ALPHA[iter] <- NB.alpha
			if (length(other.pars) > 0) OTHER.PARS[,iter] <- other.pars
		}
		if (Model=="Poisson"){
			updates <- lambda.updater(x=xx[order(x.order),],lambda=lambda[order(x.order)],lambda.additional=lambda.additional,other.pars=other.pars)
			lambda <- updates$lambda[x.order]
			lambda.additional <- updates$lambda.additional
			other.pars <- updates$other.pars
			LAMBDA[,iter] <- updates$lambda
			NB.ALPHA[iter] <- NA
			if (length(other.pars) > 0) OTHER.PARS[,iter] <- other.pars
		}
	}
	if (verbose==1) x.order <- X.ORDER
	list(X=X[,,seq(1,ndraws + burnin,by=THIN)],LAMBDA=LAMBDA[,seq(1,ncol(LAMBDA),by=THIN)],NB.ALPHA=NB.ALPHA[seq(1,length(NB.ALPHA),by=THIN)],OTHER.PARS=OTHER.PARS,x.order=x.order)
}


