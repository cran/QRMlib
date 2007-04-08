\name{GPD}
\alias{pGPD}
\alias{qGPD}
\alias{dGPD}
\alias{rGPD}
\title{
  Generalized Pareto Distribution 
}
\description{
Cumulative probability, quantiles, density and random generation from 
the generalized Pareto distribution. 
}
\usage{
pGPD(q, xi, beta=1) 
qGPD(p, xi, beta=1) 
dGPD(x, xi, beta=1, logvalue=FALSE) 
rGPD(n, xi, beta=1)
}
\arguments{
\item{x}{
vector of values at which to evaluate density 
}
\item{q}{
vector of quantiles 
}
\item{p}{
vector of probabilities 
}
\item{n}{
sample size 
}
\item{xi}{
shape parameter 
}
\item{beta}{
scale parameter 
}
\item{logvalue}{
whether or not log values of density should be returned (useful for ML) 
}
}
\value{
Probability (pGPD), quantile (qGPD), density (dGPD) or random sample (rGPD) 
for the GPD with shape xi. 
}
\seealso{
   \code{\link{fit.GPD}},
   \code{\link{pGEV}},
   \code{\link{fit.GEV}}
}
\examples{
\dontrun{
#Build a loglikelihood function for MLE which can be called from optim()
negloglik <- function(theta)
{
  -sum(dGPD(excesses.nl,theta[1],abs(theta[2]),logvalue=TRUE));
}
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

