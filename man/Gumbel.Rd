\name{Gumbel}
\alias{dGumbel}
\alias{qGumbel}
\alias{pGumbel}
\alias{rGumbel}
\title{
  Gumbel Distribution 
}
\description{
Density, quantiles, and cumulative probability of the Gumbel distribution. 
The standard Gumbel has mu value of 0 and sigma value of 1.
}
\usage{
dGumbel(x, mu=0, sigma=1, logvalue=FALSE) 
qGumbel(p, mu=0, sigma=1) 
pGumbel(q, mu=0, sigma=1)
rGumbel(n, mu=0, sigma=1)
}
\arguments{
\item{x}{
vector of values at which to evaluate density or cdf
}
\item{q}{
vector of quantiles 
}
\item{p}{
vector of probabilities 
}
\item{mu}{
location parameter of Gumbel distribution
}
\item{sigma}{
scale parameter (sigma >= 0) of Gumbel distribution
}
\item{logvalue}{
whether or not log values of density should be returned (useful for ML) 
}
\item{n}{
number of values to simulate for random Gumbel
}
}
\value{
Probability (pGumbel), quantile (qGumbel), density (dGumbel) and random vector
 (rGumbel) for the  Gumbel distribution with location parameter mu and scale
 parameter sigma. 
}
\examples{
#Simulate 1000 Gumbel variates:
rGumbelSim <- rGumbel(1000, 1.0, 2.5);
#Compare cdf of GEV to that of Gumbel when xi = 0 with location parameter 1 and scale 2.5
quantValue <- 4.5;
pG <- pGEV(q=quantValue, xi=0, mu=1.0, sigma=2.5) 
pg <- pGumbel(q=quantValue, mu=1.0, sigma=2.5);
}

\keyword{distribution}
