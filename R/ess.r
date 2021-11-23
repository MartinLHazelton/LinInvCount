#' Effective sample size
#'
#' Version of effectiveSize from coda package, modified for ease of use with DynamicLatticeBasis package.
#' @param x A matrix of MCMC out, with rows correspondings to variables and columns to iterations.
#' @return Effective sample sizes for each variable
#' @export
#' @examples
#' data(BXno3way)
#' res <- Xsampler(A=BXno3way$A,y=BXno3way$y,lambda=BXno3way$lambda,Model="Poisson",Method="Gibbs",tune.par=0.5,combine=TRUE)
#' effectiveSize(res$X)

effectiveSize <- function (x) {
  require(coda)
  if (is.mcmc.list(x)) {
    ess <- do.call("rbind", lapply(x, effectiveSize))
    ans <- apply(ess, 2, sum)
  }
  else {
    x <- as.mcmc(x)
    x <- as.matrix(x)
    variances <- apply(x, 2, var)
    indx <- variances > 0
    spec <- numeric(ncol(x))
    if (sum(indx) >0 ) spec[indx] <- spectrum0.ar(x[,indx])$spec
    ans <- ifelse(spec==0, 0, nrow(x) * variances/spec)
  }
  return(ans)
}
