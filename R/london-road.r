#' London Road network tomography data 
#'
#' Data on traffic counts observed at a sequence of points on London Road in Leicester, UK. The aim is to sample the corresponding origin-destination traffic volumes that are consistent with those counts.
#' The dataset appears as a list prepared for that purpose. 
#' Components of the list are A (the configuration matrix, which in this case is the link-path incidence matrix); y (observed traffic counts); lambda (assumed mean origin-destination traffic volumes); 
#' and MarkovBasis (a matrix with columns corresponding to the vectors in full Markov basis for the resampling problem).
#'
#' @docType data
#'
#' @usage data(LondonRoad) 
#' @examples
#' data(LondonRoad)
#' Xsampler(A=LondonRoad$A,y=LondonRoad$y,lambda=LondonRoad$lambda,Model="Poisson",Method="Gibbs",tune.par=0.5,combine=FALSE)
"LondonRoad"