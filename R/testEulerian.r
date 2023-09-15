#' Tests whether matrix is Eulerian
#'
#' Tests whether matrix is Eulerian, or c0-Eulerian, or cr0-Eulerian
#' @param U A lattice basis or matrix of signs for lattice basis
#' @param c0 Logical - test for c0-Eulerian (Eulerian with zero column sums) 
#' @param r0 Logical - test for r0-Eulerian (Eulerian with zero row sums) 
#' @param SIGN Logical - is U a matrix of signs? (If false, U is redefined as matrix of signs in function) 
#' @param verbose Monitor progress? Defaults to false.
#' @return List containing UM (indicator for unimodularity) and d (vector determinants of submatrices checked) 
#' @export

testEulerian <- function(U,c0=FALSE,r0=FALSE,SIGN=T,verbose=FALSE){
	if (SIGN) U <- sign(U)
	U <- as.matrix(U[rowSums(abs(U)) >= 2,])
	U <- as.matrix(U[,colSums(abs(U)) >= 2])
	if (c0) U <- as.matrix(U[,colSums(abs(U))!=abs(colSums(U))])
	if (r0) U <- as.matrix(U[rowSums(abs(U))!=abs(rowSums(U)),])
	Nrow <- nrow(U)
	Ncol <- ncol(U)
	if (min(Nrow,Ncol) < 2) return(FALSE)
	for (n in 2:min(Nrow,Ncol)){ 
		if (verbose) print(n)
   		colPerms<-combn(Ncol,n)
   		rowPerms<-combn(Nrow,n)
   		for (i in 1:ncol(colPerms)){
    			for (j in 1:ncol(rowPerms)){
   				subU <- U[rowPerms[,j],colPerms[,i]]
				if (any(subU!=0)){
					a <- colSums(subU)
					b <- rowSums(subU)
					if (!c0 & !r0 & all(a%%2==0) & all(b%%2==0)) return(TRUE)
					if (c0 & r0 & all(a==0) & all(b==0)) return(TRUE)
					if (c0 & !r0 & all(a==0) & all(b%%2==0)) return(TRUE)
				}
			}
		}
	}
	return(FALSE)	
}