\name{MCECMupdate}
\alias{MCECMupdate}
\title{
  MCECM Update Step for Generalized Hyperbolic 
}
\description{
updates estimates of mixing parameters in EM estimation of 
generalized hyperbolic 
}
\usage{
MCECMupdate(data, mix.pars, mu, Sigma, gamma, optpars, optfunc, xieval=FALSE)
}
\arguments{
\item{data}{
data matrix 
}
\item{mix.pars}{
current values of lambda, chi and psi 
}
\item{mu}{
current value of mu 
}
\item{Sigma}{
current value of Sigma 
}
\item{gamma}{
current value of gamma 
}
\item{optpars}{
parameters we need to optimize over (may differ from case to case) 
}
\item{optfunc}{
the function to be optimized 
}
\item{xieval}{
is it necessary to evaluate the log moment xi? 
}
}
\value{
list containing new estimates of mixing parameters as well as convergence information for optimization 
}
\details{
see Algorithm 3.14, steps (5) and (6) on page 83 of QRM 
}
\seealso{
   \code{\link{fit.mNH}}
}
\keyword{multivariate}

