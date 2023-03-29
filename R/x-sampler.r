#' Z-polytope sampler for linear inverse problems
#'
#' For the linear inverse problem y=Ax where y and x are counts, the y-fibre (solution set) is a Z-polytope (i.e. the points on the integer lattice within a convex polytope). This function implements samplers for the Z-polytope, with using a dynamic lattice basis or a full Markov basis. The underyling model for x can be Poisson, uniform or negative binomial.
#' @param y Vector of observed count data.
#' @param A Model configuration matrix, assumed to be binary.
#' @param lambda Mean vector for x.
#' @param U Optional matrix the columns of which should be a Markov (sub)-basis.
#' @param Method "MH" for Metropolis-Hastings sampler, "Gibbs" for Gibbs sampler.
#' @param Reorder Should the columns of A be reordered? Defaults to TRUE.
#' @param tune.par Tuning parameter (alpha) controlling variation in fitness values for lattice bases. Defaults to 0.5.
#' @param combine Should extra moves be included combining lattice basis vectors? Defaults to FALSE, but should usually be set to TRUE if A is not unimodular.
#' @param x.order If Reorder=FALSE, x.order can be used to reorder columns of A to match ordering of entries of x. Defaults to NULL when no such reordering is performed.
#' @param x.ini Vector of initial values for x. Default is NULL, when initial values derived through integer programming.
#' @param Model "Poisson", "Uniform", "NegBin" or "Normal" (the last being a discrete approximation).
#' @param Proposal "NonUnif" or "Unif" (default).
#' @param NB.alpha Dispersion parameter for negaqtive-binomial distribution. Defaults to 1.
#' @param ndraws Number of iterations to run sampler after burn-in. One iteration comprises cycling through the full basis (possibly augmented by a combined move). Defaults to 10^4.
#' @param burnin Number of iteractions for burn in period. Defaults to 2000, which is usually more than adequate.
#' @param verbose Controls level of detail in recording lattice bases used.
#' @param RandomMove If TRUE, sampling directions are selected at random at each iteration. If FALSE (default), the sampler cycles systematically through the available sampling directions. 
#' @param THIN Thinning parameter for output. Defaults to 1 (no thinning).
#' @return A list with components X (a matrix, each row corresponding to samplers for an entry of x) and x.order (a vector describing dynamic selection of lattice bases, if verbose=1).
#' @export
#' @examples 
#' data(LondonRoad)
#' Xsampler(A=LondonRoad$A,y=LondonRoad$y,lambda=LondonRoad$lambda,Model="Poisson",Method="Gibbs",tune.par=0.5,combine=FALSE)

Xsampler <- function (y, A, lambda, U=NULL, Method="MH", Reorder=TRUE, tune.par=0.5, combine=FALSE, x.order=NULL, x.ini=NULL, Model="Poisson", Proposal="Unif", NB.alpha=1, ndraws = 10000, burnin = 2000, verbose = 0, RandomMove=FALSE, THIN = 1) {
	require(lpSolve)
	require(numbers)
	require(extraDistr)
  	if(Model=="NegBin" & NB.alpha<=0) NB.alpha=1
	if(Model=="Uniform") Method <- "Gibbs"
	if(combine==TRUE) Proposal="Unif"
	y <- as.numeric(y)
        zero.cols.ind <- (colSums(A)==0)
        non.zero.cols.ind <- (colSums(A) > 0)
        zero.cols <- sum(zero.cols.ind)
        if (zero.cols > 0){
		r.full <- ncol(A)
		A <- A[,non.zero.cols.ind]
		lambda.full <- lambda
		lambda <- lambda[non.zero.cols.ind]
		x.ini <- x.ini[non.zero.cols.ind]
		x.order <- x.order[non.zero.cols.ind]
    	}
	r <- ncol(A)
	n <- nrow(A)
	tol <- 10^{-10}

	if (is.matrix(U)) Reorder=FALSE
	
	if (Reorder){
		lambda.star <- lambda
		lam.order <- order(lambda.star,decreasing=TRUE)  
		A <- A[,lam.order]
		x.order <- lam.order
		for (i in 3:n){
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

	X <- matrix(0, r, m*(ndraws + burnin)+1)
	if (verbose==1) X.ORDER <- matrix(0, r, m*(ndraws + burnin)+1)

	if (!is.null(x.ini)) x.ini <- x.ini[x.order]
	if (is.null(x.ini)){
		x.ini <- lp("max",objective.in=rep(1,r),const.mat=A,const.dir=rep("=",nrow(A)),const.rhs=y,all.int=T)$solution
	}
	x <- x.ini
      X[x.order,1] <- x
	if (verbose==1) X.ORDER[,1] <- x.order
	
	for (iter in seq(1, ndraws + burnin)) {
		for (jjj in 1:m){
			if (!RandomMove) j <- jjj
			if (RandomMove)  j <- sample(1:m,1)
			if (j <= ncol(U)) z <- U[,j]
			if (j > ncol(U)){
				delta <- 0
				while(sum(delta)==0) delta <- rpois(ncol(U),lambda=0.5)
				delta <- delta*sample(c(-1,1),size=ncol(U),replace=T)
				z <- U%*%delta
			}
			if (abs(dA1)!=1) { if(is_wholenumber(z)==F) z <- round(z*abs(dA1))/mGCD(round(abs(z*dA1))) }
			max.move <- floor(min((x/abs(z))[which(z<0)]))
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
			if (tune.par > 0){
                  	lambda.star <- rnorm(r,mean=lambda,sd=tune.par*lambda)
				ii <- sample(1:n,1)
				swap.indx <- abs(C[ii,])>tol
				if (any(swap.indx)){
					jj <- sample((1:(r-n))[swap.indx],1)
					if (lambda.star[ii] <= lambda.star[jj+n]/abs(C[ii,jj])){
					ei <- rep(0,n)
 					ej <- rep(0,r-n)
					ei[ii] <- 1
					ej[jj] <- 1
					dA1 <- round(C[ii,jj]*dA1)
					C <- C - outer(C[,jj]-ei,C[ii,]+ej)/C[ii,jj]
					U <- rbind(-C,diag(r-n))
					x.order[c(ii,jj+n)] <- x.order[c(jj+n,ii)]
					x[c(ii,jj+n)] <- x[c(jj+n,ii)]
 					lambda[c(ii,jj+n)] <- lambda[c(jj+n,ii)]
 					}
				}
			}
			X[x.order,1+(iter-1)*m+jjj] <- x
			if (verbose==1) X.ORDER[,1+(iter-1)*m+jjj] <- x.order
		}
	}
	if (verbose==1) x.order <- X.ORDER
    	if (zero.cols > 0){
		XX <- X
		X <- matrix(NA,nrow=r.full,ncol=ncol(XX))
		X[non.zero.cols.ind,] <- XX
		for (i in 1:zero.cols){
			X[which(zero.cols.ind)[i],] <- rpois(ncol(XX),lambda.full[which(zero.cols.ind)[i]])
		}
   	}
	list(X=X[,seq(1,ncol(X),by=THIN)],x.order=x.order)
}


