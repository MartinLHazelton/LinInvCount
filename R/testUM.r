#' Tests unimodularity for configutation matrix
#'
#' Tests unimodularity for configutation matrix. For large matrices, checks just a sample of submatrices.
#' @param A A configuration matrix
#' @param sim.size Just check a subset of sim.size submatrices if number of submatrices exceeds that number. Defaults to 1e6.
#' @param verbose Monitor progress? Defaults to false.
#' @return List containing UM (indicator for unimodularity) and d (vector determinants of submatrices checked) 
#' @export

testUM <- function(A,sim.size=1e6,verbose=FALSE){
	n <- nrow(A)
	r <- ncol(A)
	if (choose(r,n) <= sim.size) colPerms <- combn(r,n)
	if (choose(r,n) > sim.size) colPerms <- replicate(sim.size,sample(1:r,size=n,replace=F))
	dd <- numeric(ncol(colPerms))
	UM <- TRUE
	for (i in 1:ncol(colPerms)){
		if (verbose & i%%100==0) cat(i,"\n")
		d <- det(A[,colPerms[,i]])
		if (abs(d) > 1) UM <- FALSE
		dd[i] <- d
	}
	list(UM=UM,d=dd)	
}