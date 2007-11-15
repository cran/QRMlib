\name{fit.binomial}
\alias{fit.binomial}
\title{
  Fit Binomial Distribution 
}
\description{
fits binomial distribution by maximum likelihood 
}
\usage{
fit.binomial(M, m)
}
\arguments{
\item{M}{
vector of numbers of successes 
}
\item{m}{
vector of numbers of trials 
}
}
\value{
list containing parameter estimates and details of fit 
}
\seealso{
   \code{\link{fit.binomialBeta}},
   \code{\link{fit.binomialLogitnorm}},
   \code{\link{fit.binomialProbitnorm}}
}

\examples{
data(spdata.raw); 
attach(spdata.raw); 
BdefaultRate <- Bdefaults/Bobligors; 
mod0 <- fit.binomial(Bdefaults, Bobligors); 
}
\keyword{methods}

