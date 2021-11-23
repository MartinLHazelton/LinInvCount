#' Yang's network tomography data 
#'
#' Simulated traffic counts on Yang's (1995) network. The aim is to sample the corresponding origin-destination traffic volumes that are consistent with those counts.
#' The dataset appears as a list prepared for that purpose. 
#' Components of the list are A (the configuration matrix, which in this case is the link-path incidence matrix); y (observed traffic counts); lambda (assumed mean origin-destination traffic volumes); 
#' and MarkovBasis (a matrix with columns corresponding to the vectors in full Markov basis for the resampling problem).
#'
#' @references Yang, H. (1995). Heuristic algorithms for the bilevel origin-destination matrix estimation problem. Transportation Research Part B: Methodological, 29(4), 231-242.
#'
#' @docType data
#'
#' @usage data(YangNetwork) 
#' @examples
#' data(YangNetwork)
#' Xsampler(A=YangNetwork$A,y=YangNetwork$y,lambda=YangNetwork$lambda,Model="Poisson",Method="Gibbs",tune.par=0.5,combine=TRUE)
"YangNetwork"