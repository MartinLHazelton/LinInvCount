#' Remove redundancy frrom linear system
#'
#' Removes redundant rows of configuration matric and corresponding entries in observed data
#' @param A Configuration matrix
#' @param y Vector of observed data. Should be of length nrow(A)
#' @return A list with components A and y, respectively the configuration matrix and data vector with redundant rows removed.
#' @export

remove_redundancy <- function(A,y){
	if (length(y)!=nrow(A)) stop("A and y incompatible")
	tmp <- qr(t(A))
	n <-tmp$rank
	indx <- tmp$pivot[1:n]
	A <- A[indx,]
	y <- y[indx]
	list(A=A,y=y)
}