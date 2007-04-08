\name{claytonmix}
\alias{claytonmix}
\alias{dclaytonmix}
\alias{pclaytonmix}
\alias{rclaytonmix}
\title{
  Mixing Distribution on Unit Interval Yielding Clayton Copula Model 
}
\description{
density, cumulative probability, and random generation for a mixture 
distribution on the unit interval which gives an exchangeable Bernoulli 
mixture model equivalent to a Clayton copula model 
}
\usage{
dclaytonmix(x, pi, theta) 
pclaytonmix(q, pi, theta) 
rclaytonmix(n, pi, theta)
}
\arguments{
\item{x}{
values at which density should be evaluated 
}
\item{q}{
values at which cumulative distribution should be evaluated 
}
\item{n}{
sample size 
}
\item{pi}{
parameter of distribution 
}
\item{theta}{
parameter of distribution 
}
}
\value{
values of density (dclaytonmix), distribution function (pclaytonmix) 
or random sample (rclaytonmix) 
}
\details{
see page 362 in QRM 
}
\seealso{
   \code{\link{dbeta}},
   \code{\link{dprobitnorm}}
}
\examples{
#probability of only one obligor defaulting B class (see Table 8.6 in QRM book)
pi.B <- 0.0489603; 
#joint probability of two obligors defaulting B class (see Table 8.6 in QRM book)
pi2.B <- 0.003126529; 
# Calibrate Calyton copula model to pi.B and pi2.B
claytonmix.pars <- cal.claytonmix(pi.B,pi2.B)
# We could also look at mixing densities. Get probability of Clayton mix
# This picture essentially shows large sample asymptotics
#Build 1000 equally-spaced values on unit interval (multiples of .000999); 
#discard all values except those below 0.25
q <- (1:1000)/1001;
q <- q[q<0.25]; #reduce to lowest 250 values
#get probabilities for each of 250 lowest values on unit interval
d.claytonmix <- dclaytonmix(q,claytonmix.pars[1],claytonmix.pars[2]);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}
