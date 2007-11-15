\name{GEV}
\alias{pGEV}
\alias{qGEV}
\alias{dGEV}
\alias{rGEV}
\title{
  Generalized Extreme Value Distribution 
}
\description{
Cumulative probability, quantiles, density and random generation from 
the generalized extreme value distribution. 
}
\usage{
pGEV(q, xi, mu=0, sigma=1) 
qGEV(p, xi, mu=0, sigma=1) 
dGEV(x, xi, mu=0, sigma=1, logvalue=FALSE) 
rGEV(n, xi, mu=0, sigma=1)
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
\item{mu}{
location parameter 
}
\item{sigma}{
scale parameter 
}
\item{logvalue}{
whether or not log values of density should be returned (useful for ML) 
}
}
\value{
Probability (pGEV), quantile (qGEV), density (dGEV) or random sample (rGEV) 
for the GEV distribution with shape xi (with location parameter mu and 
location parameter sigma)
}
\seealso{
   \code{\link{fit.GEV}},
   \code{\link{fit.GPD}},
   \code{\link{pGPD}}
}
\examples{
#Compare cdf of GEV to that of Gumbel when xi = 0 with location parameter 1 and scale 2.5
quantValue <- 4.5;
pG <- pGEV(q=quantValue, xi=0, mu=1.0, sigma=2.5) 
pg <- pGumbel(q=quantValue, mu=1.0, sigma=2.5);
}
\keyword{distribution}

