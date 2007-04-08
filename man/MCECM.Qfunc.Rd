\name{MCECM.Qfunc}
\alias{MCECM.Qfunc}
\title{
  Optimization Function for MCECM Fitting of GH 
}
\description{
a functional form that must be optimized when fitting members of 
generalized hyperbolic family with an MCECM algorithm 
}
\usage{
MCECM.Qfunc(lambda, chi, psi, delta, eta, xi)
}
\arguments{
\item{lambda}{
lambda parameter 
}
\item{chi}{
chi parameter 
}
\item{psi}{
pi parameter 
}
\item{delta}{
delta statistic 
}
\item{eta}{
eta statistic 
}
\item{xi}{
xi statistic 
}
}
\value{
value of function 
}
\details{
this is the Q2 function on page 82 of QRM 
}
\seealso{
   \code{\link{MCECMupdate}},
   \code{\link{EMupdate}},
   \code{\link{fit.mNH}}
}
\keyword{multivariate}

