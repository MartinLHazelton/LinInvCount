#' Maximum edge loading
#'
#' Computes maximum edge loading for transition graph.
#' @param A A configuration matrix
#' @param U A paritition lattice basis (given by the columns of the matrix)
#' @param y Vector of observed counts
#' @param U.reorder Recordering of columns of A/rows of U
#' @param lambda Rate parameters, need specifying if using Poisson model
#' @param Fibre Fibre (solution set) for y=Ax
#' @param UNIFORM Indicator for uniform model (otherwise Poisson)
#' @return List containing loading and identifier for most loaded edge
#' @export

MaxEdgeLoading <- function(A,U,y,U.reorder=1:nrow(U),lambda=NULL,Fibre=NULL,UNIFORM=TRUE){
	POISSON <- !UNIFORM
	if (POISSON & is.null(lambda)) stop("Non-null lambda needed for Poisson models")
	U <- U[U.reorder,]
	if (is.null(Fibre)) Fibre <- FindFibre(A,y)
	r <- ncol(A)
	n <- nrow(A)
	N <- ncol(Fibre)
	transition.matrix <- matrix(0,N,N)
	d <- ncol(U)
	for (from in 1:N){
		for (u in 1:d){
			x <- Fibre[,from]
			z <- U[,u]
           		max.move <- floor(min((x/abs(z))[which(z < 0)]))
            	min.move <- -floor(min((x/abs(z))[which(z > 0)]))
			ray.length <- max.move-min.move+1
			if (ray.length==1) transition.matrix[from,from] <- transition.matrix[from,from]+1
			if (ray.length > 1){
           			x.min <- x + min.move * z
           			x.max <- x + max.move * z
         			if (POISSON){
					x.matrix <- round(t(mapply(seq, from = x.min, by = z, length.out = ray.length)))
 					probs <- exp(colSums(dpois(x.matrix, lambda, log = T)))
					probs <- probs/sum(probs)
				}
				for (step in 1:ray.length){
					xto <- x.min + z*(step-1)
					to <- which(colSums(abs(Fibre-xto))==0)
					if(UNIFORM) transition.matrix[from,to] <- transition.matrix[from,to] + 1/ray.length
					if(POISSON) transition.matrix[from,to] <- transition.matrix[from,to] + probs[step]	
				}	
			}
		}
	}
	transition.matrix <- transition.matrix/d
	canonical.paths <- matrix(0,nrow=d,ncol=N^2)
	prob.product <- numeric(N^2)
	if (POISSON) normalizing.constant <- sum(apply(dpois(Fibre,lambda),2,prod)) 
	for (from in 1:N){
		for (to in 1:N){
			z <- Fibre[,to] - Fibre[,from]
			canonical.paths[,(from-1)*N+to] <- z[-(1:n)]
			if (UNIFORM) prob.product[(from-1)*N+to] <- 1/N^2
			if (POISSON) prob.product[(from-1)*N+to] <- prod(dpois(Fibre[,from],lambda))*prod(dpois(Fibre[,to],lambda))/normalizing.constant^2
		}
	}
	rho <- 0
	for (from in 1:N){
		for (to in 1:N){
			if(transition.matrix[from,to]>0 & from !=to){
				edge <- Fibre[,to] - Fibre[,from]
				edge2 <- edge[-(1:n)]
				indx <- colSums(1*((canonical.paths[edge2!=0,,drop=FALSE]/edge2[edge2!=0])==1))==sum(edge2!=0)
				#cat(prob.product[indx],transition.matrix[from,to],"\n")
				edge.loading <- sum(prob.product[indx])/transition.matrix[from,to]
				if (edge.loading > rho){
					rho <- edge.loading
					max.edge <- edge
				}
			}
		}
	}
	list(MaxEdgeLoading=rho,MostLoadedEdge=max.edge)
}
