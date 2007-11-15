\name{probitnorm}
\alias{probitnorm}
\alias{dprobitnorm}
\alias{pprobitnorm}
\alias{rprobitnorm}
\title{
  Probit-Normal Distribution 
}
\description{
density, cumulative probability and random number generation for distribution of random variable Q on unit interval such that the probit transform of Q has a normal distribution with parameters mu and sigma 
}
\usage{
dprobitnorm(x, mu, sigma) 
pprobitnorm(q, mu, sigma) 
rprobitnorm(n, mu, sigma)
}
\arguments{
\item{x}{
vector of values in unit interval at which to evalualte density 
}
\item{q}{
vector of values in unit interval at which to evalualte cumulative probabilities 
}
\item{n}{
sample size 
}
\item{mu}{
scalar parameter 
}
\item{sigma}{
scalar parameter 
}
}
\value{
vector of density values (dprobitnorm), cumulative probabilities (pprobitnorm) or random sample (rprobitnorm) 
}
\details{
see pages 353-354 in QRM 
}
\seealso{
   \code{\link{dbeta}},
   \code{\link{dclaytonmix}}
}
\examples{
#MODEL RISK See especially Section 8.4.6 on p. 364 of QRM book
data(spdata.raw);
attach(spdata.raw);
pi.B <- momest(Bdefaults, Bobligors)[1];#one obligor defaulting pi = .04896
#second value returned is probability of joint default probability for two firms.
pi2.B <- momest(Bdefaults, Bobligors)[2]; #two obligors defaulting jointly pi2 = .0031265
#Build 1000 equally-spaced value on unit interval as multiples of .000999; discard 
#all values except those below 0.25 because we want to look at the tail, i.e. Q > 0.25 
#via the tail function [1 - P(Q <= 0.25)]
# Model Risk Experiment
# Calibrate a 1-Factor Creditmetrics (probitnormal) model to pi.B and pi2.B for all models:
#The following values are shown in Table 8.6, column B, row labeled 'Probit-normal'.  
#In other words, find the probitnorm mu and sigma values which give same probabilities as 
#momest()
probitnorm.pars <- cal.probitnorm(pi.B,pi2.B);
probitnorm.pars;
q <- (1:1000)/1001;
q <- q[q<0.25];
# We could also look at mixing densities. Remember that density values for continuous 
#variables may exceed 1 since they give an approximation for the change in the cdf value 
#as we change the x value. Hence if the cdf increases by 0.2 as we increase x from 0.1 to 
#0.2, the density should be about 2.0 (dF(x)/dx).
d.probitnorm <- dprobitnorm(q,probitnorm.pars[1],probitnorm.pars[2]);
}

\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

