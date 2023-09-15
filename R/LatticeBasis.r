#' Computes partition lattice basis
#'
#' Computes solution set (fibre) for y=Ax for non-negative integers
#' @param A A configuration matrix
#' @param mu Optional specification of E[x] promote favourable ordering of colums of A (i.e. entries of x)
#' @param reorder All columns of A to be reordered to compute basis? Output still matches original column ordering
#' @return List containing partition lattice basis U, determinant of A1 in A = [A1, A2], column ordering
#' @export

LatticeBasis <- function(A,mu=NULL,reorder=T){
      r <- ncol(A)
	n <- nrow(A)
	if (is.null(mu)) mu <- rnorm(r)
      col.order <- order(mu, decreasing = TRUE)
      if (!reorder) col.order <- 1:r 
	A <- A[, col.order]
      for (i in 2:n) {
          while (qr(A[, 1:i])$rank < i) {
              A <- A[, c(1:(i - 1), (i + 1):r, i)]
              col.order <- col.order[c(1:(i - 1), (i + 1):r, i)]
          }
      }
	A1 <- A[, 1:n]
      dA1 <- det(A1)
      A2 <- as.matrix(A[, -c(1:n)])
      A1.inv <- solve(A1)
      C <- A1.inv %*% A2
      U <- rbind(-C, diag(r - n))
	U[col.order,] <- U
	list(U=U,detA1=dA1,col.order=col.order)	
}