\name{ghypB}
\alias{dghypB}
\alias{rghypB}
\title{
  Univariate Generalized Hyperbolic Distribution B 
}
\description{
Density and random number generation for univariate generalized hyperbolic 
distribution in standard parameterization (alpha-beta-delta). (The dispersion 
matrix Sigma is identically 1, i.e. a scalar 1.) See pp. 77-81 in QRM. 
}
\details{
See page 78 in QRM for joint density formula (3.30) with Sigma (dispersion matrix)
the identity and d=1 (meaning a univariate distribution) applies. 

The B parameterization corresponds to the original alpha-beta-delta model used by 
Blaesild (1981) in earlier literature.  If gamma is 0, we have a normal variance 
mixture defined by the paramters alpha-beta-delta.  This thickens the tail.

If gamma exceeds zero, we have a normal mean-variance mixture where the mean 
is also perturbed to equal mu + (W * gamma) which introduces ASYMMETRY as well. 

Values for lambda and mu are identical in both QRM and B parameterizations.

Sigma does not appear in parameter list since in the univariate case its value is 
assumed to be identically 1.
}
\note{
Density values from dgyhp() should be identical to those from dghypB() 
if the alpha-beta-delta parameters of the B type are translated to the
corresponding gamma-chi-psi parameters of the QRM type by formulas on pp 79-80.
}

\usage{
dghypB(x, lambda, delta, alpha, beta=0, mu=0, logvalue=FALSE) 
rghypB(n, lambda, delta, alpha, beta=0, mu=0)
}
\arguments{
\item{x}{
values at which to evaluate density 
}
\item{n}{
sample size 
}
\item{lambda}{
scalar parameter 
}
\item{delta}{
scalar parameter 
}
\item{alpha}{
scalar parameter 
}
\item{beta}{
skewness parameter 
}
\item{mu}{
location parameter 
}
\item{logvalue}{
Should log density be returned? Default is FALSE 
}
}
\value{
values of density or log-density (dghypB) or random sample (rghypB) 
}
\seealso{
   \code{\link{dghyp}},
   \code{\link{besselM3}}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

