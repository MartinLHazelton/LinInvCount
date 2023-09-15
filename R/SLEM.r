#' SLEM
#'
#' Computes second largest eigvenvalue modulus for transition matrix
#' @param A A configuration matrix
#' @param U A paritition lattice basis (given by the columns of the matrix)
#' @param y Vector of observed counts
#' @param lambda Rate parameters, need specifying if using Poisson model
#' @param Fibre Fibre (solution set) for y=Ax
#' @param UNIFORM Indicator for uniform model (otherwise Poisson)
#' @return List containing SLEM, transition matrix and fibre
#' @export

SLEM <- function(A,U,y,lambda=NULL,Fibre=NULL,UNIFORM=TRUE){
	POISSON <- !UNIFORM
	if (POISSON & is.null(lambda)) stop("Non-null lambda needed for Poisson models")
	if (is.null(Fibre)) Fibre <- FindFibre(A,y)
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
	slem <- max(Mod(eigen(transition.matrix)$values[2]),Mod(eigen(transition.matrix)$values[length(eigen(transition.matrix)$values)]))
	list(SLEM=slem,transition.matrix=transition.matrix,Fibre=Fibre)
}
