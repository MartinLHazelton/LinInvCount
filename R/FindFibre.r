#' Find fibre
#'
#' Computes solution set (fibre) for y=Ax for non-negative integers
#' @param A A configuration matrix
#' @param y Vector of observed counts
#' @return Matrix with columns giving fibre elements
#' @export

FindFibre <- function(A,y){
      r <- ncol(A)
	n <- nrow(A)
	col.order <- 1:r 
      for (i in 2:n) {
          while (qr(A[, 1:i])$rank < i) {
              A <- A[, c(1:(i - 1), (i + 1):r, i)]
              col.order <- col.order[c(1:(i - 1), (i + 1):r, i)]
          }
      }
	Apart1 <- A[,1:n]
	Apart1inv <- solve(Apart1)
	Apart2 <- A[,-(1:n)]
	bounds <- y*Apart2 + 1/Apart2 - Apart2
	xmax <- apply(bounds,2,min)
	Fibre <- numeric(0)
	gridlist <- vector(mode = "list", length = r-n)
	for (i in 1:(r-n)) gridlist[[i]] <- 0:xmax[i]
	grid <- expand.grid(gridlist)
	for (i in 1:nrow(grid)){
		x2 <- as.numeric(grid[i,])
		x1 <- Apart1inv%*%(y-Apart2%*%x2)
		x <- c(x1,x2)
		if (all(x>=0)) Fibre <- cbind(Fibre,x)
	}
	Fibre[col.order,] <- Fibre
	Fibre
}