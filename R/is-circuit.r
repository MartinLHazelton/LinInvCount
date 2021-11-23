#' Check whether basis vectors are circuits 
#'
#' A circuit is a vector with minimal support by inclusion and coprime (non-zero) entries. This function checks whether each column of a matrix U is a circuit of the solution set for y = Ax
#' @param A Configuration matrix with ncol(A) >= nrow(A)
#' @param U A matrix with columns comprising some basis for the integer lattice kernel of A
#' @return A logical vector of length ncol(U) indicating whether each column of U is a circuit
#' @export

is_circuit <- function(A,U){
	r <- nrow(U)
	nU <- ncol(U)
	circuit <- numeric(nU)
	for (i in 1:nU){
		supp <- (U[,i] != 0)
		Ai <- as.matrix(A[,supp])
		circuit[i] <- (1+qr(Ai)$rank)==nrow(Ai)
	}
	circuit
}
