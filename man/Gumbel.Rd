\name{Gumbel}
\alias{dGumbel}
\alias{qGumbel}
\alias{pGumbel}
\title{
  Gumbel Distribution 
}
\description{
Density, quantiles, and cumulative probability of the Gumbel distribution. 
}
\usage{
dGumbel(x, logvalue=FALSE) 
qGumbel(p) 
pGumbel(q)
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
\item{logvalue}{
whether or not log values of density should be returned (useful for ML) 
}
}
\value{
Probability (pGumbel), quantile (qGumbel) and density (dGumbel) for the 
Gumbel distribution with shape xi. 
}
\examples{
\dontrun{
#If the shape parameter of the GEV is xi=0, then the distribution collapses to
#the Gumbel with parameters properly scaled:
x <- 2.5; mu <- 1.0; sigma = 2.5;
xx = x - mu/sigma;
if(xi == 0) dGEVvalue <- dGumbel(xx, logvalue=TRUE) - log(sigma);
}
}

\keyword{distribution}
