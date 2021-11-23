#' Book-crossing 2-way contingency table data 
#'
#' 30x15 contingency table classifying book reviews by book ISBN (30 levels) and country of reviewer (15 levels).
#' The dataset appears as a list prepared for contingency table resampling conditional on the marginal totals. 
#' Components of the list are A (the configuration matrix); y (marginal totals); x (true cell counts); lambda (cell means under independence model); and MarkovBasis (a matrix with columns corresponding
#' to the vectors in full Markov basis for the resampling problem).
#'
#' @docType data
#'
#' @usage data(BX2way) 
#' @examples
#' data(BX2way)
#' Xsampler(A=BX2way$A,y=BX2way$y,lambda=BX2way$lambda,Model="Uniform",Method="Gibbs",tune.par=0.5,combine=FALSE,burnin=50,ndraws=200)
"BX2way"