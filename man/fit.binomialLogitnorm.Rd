\name{fit.binomialLogitnorm}
\alias{fit.binomialLogitnorm}
\title{
  Fit Logitnormal-Binomial Distribution 
}
\description{
fits a mixed binomial distribution where success probability has a logitnormal distribution 
This function has been altered in the R-language edition to contain two extra parameters
providing upper and lower limits for the input parameters M and m.  if convergence occurs at 
an endpoint of either limit, you need to reset lower and upper parameter estimators and run 
the function again
}
\usage{
fit.binomialLogitnorm(M, m, startvals=c(-1, 0.5), 
lowerParamLimits = c(-5.0, 0.02), upperParamLimits = c(1,0.9))
}
\arguments{
\item{M}{
vector of numbers of successes (e.g. number of defaults in a credit-rating class)
}
\item{m}{
vector of numbers of trials (e.g. number of obligors in a credit-rating class)
}
\item{startvals}{
starting values for parameter estimates 
}
\item{lowerParamLimits}{
vector with lower limits for each parameter to be used by optimization algorithm
}
\item{upperParamLimits}{
vector with upper limits for each parameter to be used by optimization algorithm
}
}
\details{
This function calls the R-language method optim(...method="L-BFGS-B") which uses input 
parameter vectors of upper and lower limits.  Hence if convergence occurs at an endpoint 
of either limit, you may need to expand the corresponding upper or lower limit and run 
the function again.
}
\value{
list containing parameter estimates and details of fit: 
\item{par.ests }{vector of optimum parameter estimators}
\item{maxloglik }{value of likelihood function at optimum}
\item{converged }{T or F indicating convergence}
\item{details }{any messages associated with convergence algorithm}
\item{pi }{probability of a single default (see p. 345 in QRM)}
\item{pi2 }{probability of two joint defaults (see p. 345 in QRM)}
\item{rhoY}{default correlation (see p. 345 in QRM)}
}
\seealso{
   \code{\link{fit.binomial}},
   \code{\link{fit.binomialBeta}},
   \code{\link{fit.binomialProbitnorm}}
}
\examples{
data(spdata.raw); 
attach(spdata.raw); 
BdefaultRate <- Bdefaults/Bobligors; 
# A little patience is required for the next model ... 
mod3 <- fit.binomialLogitnorm(Bdefaults, Bobligors); 
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

