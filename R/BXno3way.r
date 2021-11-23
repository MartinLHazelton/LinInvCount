#' Book-crossing 3-way contingency table data 
#'
#' 3x4x6 contingency table classifying book reviews by age of reviewer (three categories, 0-25, 26-50, 51+), country of reviewer (Australia, Canada, UK, US) and book (labelled A-F).
#' The dataset appears as a list prepared for contingency table resampling conditional on the sufficient statistics for a log-linear model with all 2-way but no 3-way interactions). 
#' Components of the list are A (the configuration matrix); y (sufficient statistics); x (true cell counts); lambda (cell means); and MarkovBasis (a matrix with columns corresponding
#' to the vectors in full Markov basis for the resampling problem).
#'
#' @docType data
#'
#' @usage data(BXno3way) 
#' @examples
#' data(BXno3way)
#' Xsampler(A=BXno3way$A,y=BXno3way$y,lambda=BXno3way$lambda,Model="Poisson",Method="Gibbs",tune.par=0.5,combine=TRUE)
"BXno3way"