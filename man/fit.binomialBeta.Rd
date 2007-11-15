\name{fit.binomialBeta}
\alias{fit.binomialBeta}
\title{
  Fit Beta-Binomial Distribution to defaults and obligors
}
\description{
fit a beta-binomial distribution by maximum likelihood 
}
\usage{
fit.binomialBeta(M, m, startvals=c(2, 2), ses=FALSE)
}
\arguments{
\item{M}{
vector of numbers of successes (defaults)
}
\item{m}{
vector of numbers of trials (obligors). M and m vectors will have equal length
which represents the number of different credit classifications/ratings
}
\item{startvals}{
starting values for parameter estimates 
}
\item{ses}{
whether standard errors should be calculated 
}
}
\value{
list containing parameter estimates and details of fit 
}
\seealso{
   \code{\link{fit.binomial}},
   \code{\link{fit.binomialLogitnorm}},
   \code{\link{fit.binomialProbitnorm}}
}
\examples{
data(spdata.raw);
spdata.raw; 
#attach data so we don't have to qualify the data column names:
attach(spdata.raw); 
BdefaultRate <- Bdefaults/Bobligors; 
mod1 <- fit.binomialBeta(Bdefaults, Bobligors); 
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

