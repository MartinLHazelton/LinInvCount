#' Front ends for Xsampler 
#'
#' Front ends for Xsampler, for quick implementation of dynamic lattice basis and full Markov basis samplers. XGibbs and XMH are Gibbs and Metropolis-Hastings versions of Xsampler using dynamic lattice bases. XMarkovBasis_Gibbs and XMarkovBasis_MH are versions using a full Markov basis. 
#' @param ... As for Xsampler
#' @return A list as for Xsampler
#' @export
#' @examples
#' data(BXno3way)
#' XGibbs(BXno3way$y,BXno3way$A,BXno3way$lambda,combine=T,tune.par=0.5)

XGibbs <- function (y, A, lambda, Reorder=TRUE, tune.par=0.5, combine=F, x.order=NULL, x.ini=NULL, Model="Poisson", NB.alpha=0, ndraws = 10000, burnin = 2000, verbose = 0, THIN = 1) {
	Xsampler(y=y,A=A,lambda=lambda,U=NULL,Method="Gibbs",combine=combine,Reorder=Reorder,tune.par=tune.par,x.order=x.order,x.ini=x.ini,Model=Model,Proposal=NULL,NB.alpha=NB.alpha,ndraws=ndraws,burnin=burnin,verbose=verbose,THIN=THIN) 
}


#' Front ends for Xsampler 
#'
#' Front ends for Xsampler, for quick implementation of dynamic lattice basis and full Markov basis samplers. XGibbs and XMH are Gibbs and Metropolis-Hastings versions of Xsampler using dynamic lattice bases. XMarkovBasis_Gibbs and XMarkovBasis_MH are versions using a full Markov basis. 
#' @param ... As for Xsampler
#' @return A list as for Xsampler
#' @export
#' @examples
#' data(BXno3way)
#' XGibbs(BXno3way$y,BXno3way$A,BXno3way$lambda,combine=T,tune.par=0.5)

XMH <- function (y, A, lambda, Reorder=TRUE, tune.par=0.5, combine=F, Proposal="NonUnif",x.order=NULL, x.ini=NULL, Model="Poisson", NB.alpha=0, ndraws = 10000, burnin = 2000, verbose = 0, THIN = 1) {
	Xsampler(y=y,A=A,lambda=lambda,U=NULL,Method="MH",combine=combine,Reorder=Reorder,tune.par=tune.par,x.order=x.order,x.ini=x.ini,Model=Model,Proposal=Proposal,NB.alpha=NB.alpha,ndraws=ndraws,burnin=burnin,verbose=verbose,THIN=THIN) 
}


#' Front ends for Xsampler 
#'
#' Front ends for Xsampler, for quick implementation of dynamic lattice basis and full Markov basis samplers. XGibbs and XMH are Gibbs and Metropolis-Hastings versions of Xsampler using dynamic lattice bases. XMarkovBasis_Gibbs and XMarkovBasis_MH are versions using a full Markov basis. 
#' @param ... As for Xsampler
#' @return A list as for Xsampler
#' @export
#' @examples
#' data(BXno3way)
#' XGibbs(BXno3way$y,BXno3way$A,BXno3way$lambda,combine=T,tune.par=0.5)

XMarkovBasis_Gibbs <- function (y, A, lambda, U, x.ini=NULL, Model="Poisson", NB.alpha=0, ndraws = 10000, burnin = 2000, verbose = 0, THIN = 1) {
	Xsampler(y=y,A=A,lambda=lambda,U=U,Method="Gibbs",Reorder=FALSE,x.order=NULL,x.ini=x.ini,Model=Model,Proposal=NULL,NB.alpha=NB.alpha,ndraws=ndraws,burnin=burnin,verbose=verbose,THIN=THIN) 
}


#' Front ends for Xsampler 
#'
#' Front ends for Xsampler, for quick implementation of dynamic lattice basis and full Markov basis samplers. XGibbs and XMH are Gibbs and Metropolis-Hastings versions of Xsampler using dynamic lattice bases. XMarkovBasis_Gibbs and XMarkovBasis_MH are versions using a full Markov basis. 
#' @param ... As for Xsampler
#' @return A list as for Xsampler
#' @export
#' @examples
#' data(BXno3way)
#' XGibbs(BXno3way$y,BXno3way$A,BXno3way$lambda,combine=T,tune.par=0.5)

XMarkovBasis_MH <- function (y, A, lambda, U, x.ini=NULL, Model="Poisson", Proposal="NonUnif", NB.alpha=0, ndraws = 10000, burnin = 2000, verbose = 0, THIN = 1) {
	Xsampler(y=y,A=A,lambda=lambda,U=U,Method="MH",Reorder=FALSE,x.order=NULL,x.ini=x.ini,Model=Model,Proposal="Unif",NB.alpha=NB.alpha,ndraws=ndraws,burnin=burnin,verbose=verbose,THIN=THIN) 
}