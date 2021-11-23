#' Check whether partitions of the figuration matrix are unimodular 
#'
#' Check whether partitions of the figuration matrix are unimodular. For partitions A=[A1 A2], the function um_partition computes the determinants of A1 for all possible (unique) combinations of rows. Function um_partition_sample coputes determinants of A1 for a sample of those partitions.
#' @param A Configuration matrix with ncol(A) >= nrow(A)
#' @return A list comprising a table of determinants and the number of rows (n) and columns (r) of the matrix A
#' @export
#' @examples
#' data(YangNetwork)
#' set.seed(2020)
#' um_partition_sample(YangNetwork$A)

um_partition <- function(A){
	require(gtools)
	n <- qr(t(A))$rank
	A <- t( t(A)[,qr(t(A))$pivot[1:n]] )
	r <- ncol(A)
	combs <- combinations(r,n,1:r)
	dets <- numeric(nrow(combs))
	for (i in 1:nrow(combs)){
		if (i/100000 == round(i/100000)) cat(i,"\n")
		dets[i] <- det(A[,combs[i,1:n]])
	}
	list(dets=table(dets),r=r,n=n)
}

um_partition_sample <- function(A,nsample=1e6){
	r <- ncol(A)
	n <- nrow(A)
	dets <- numeric(nsample)
	for (i in 1:nsample){
		sample.rows <- sample(1:r,n,replace=F)
		dets[i] <- det(A[,sample.rows])
	}
	list(dets=table(dets),r=r,n=n)
}