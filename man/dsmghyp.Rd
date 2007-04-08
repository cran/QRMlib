\name{dsmghyp}
\alias{dsmghyp}
\title{
  Symmetric Multivariate Generalized Hyperbolic Distribution 
}
\description{
Density of elliptical subfamily of multivariate generalized hyperbolic family.
The symmetric family is a normal-variance mixture since the gamma parameter
associated with the mean mixture is by assumption equal to zero.
}
\details{
See page 78 in QRM for joint density formula (3.30) with Sigma a d-dimensional
dispersion matrix (d > 1) consistent with a multivariate distribution). This is 
a more intuitive parameterization of the alpha-beta-delta model used by Blaesild 
(1981) in earlier literature since it associates all parameters with mixtures of
both mean and variance.  Since gamma is 0, we have a normal-variance mixture where 
the mixing variable W has a GIG generalized inverse gaussian) distribution with 
parameters lambda, chi, psi. This thickens the tail.

Since gamma equals zero, we have no perturbation of the mean so no ASYMMETRY 
is introduced and hence the distribution is symmetric.

There is no random number generation associated with the multivariate model in this 
implementation of the R-language and S-Plus code. 

See pp. 77-81 of QRM and appendix A.2.5 for details.

dsmghyp() is frequently called from the function dmghyp().
}

\usage{
dsmghyp(x, lambda, chi, psi, mu, Sigma, logvalue=FALSE)
}
\arguments{
\item{x}{
matrix with n rows and d columns; density is evaluated at each vector of row values 
}
\item{lambda}{
scalar parameter 
}
\item{chi}{
scalar parameter 
}
\item{psi}{
scalar parameter 
}
\item{mu}{
location vector 
}
\item{Sigma}{
dispersion matrix 
}
\item{logvalue}{
should log density be returned; default is FALSE 
}
}
\value{
vector of length n containing values of density or log-density 
}
\seealso{
   \code{\link{dmghyp}}
}
\examples{
\dontrun{
dmghyp <- function(x, lambda, chi, psi, mu, Sigma, gamma,logvalue=FALSE)
{
  #Call symmetric form if gamma vector is identically zero:
  if (sum(abs(gamma))==0) 
       out <- dsmghyp(x, lambda,chi,psi,mu,Sigma,logvalue=TRUE);
  # lines removed here
}
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

