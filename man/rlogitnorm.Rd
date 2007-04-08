\name{rlogitnorm}
\alias{rlogitnorm}
\title{
  Random Number Generation from Logit-Normal Distribution 
}
\description{
Random number generation for distribution of random variable Q on unit interval 
such that the probit transform of Q has a normal distribution with parameters 
mu and sigma 
}
\usage{
rlogitnorm(n, mu, sigma)
}
\arguments{
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
random sample of size n 
}
\details{
see pages 353-354 in QRM 
}
\seealso{
   \code{\link{rbeta}},
   \code{\link{rclaytonmix}},
   \code{\link{rprobitnorm}}
}
\examples{
#set number, mean, and variance:
num <- 1000;
mu <- 2.0;
sigma <- 1.25;
#Simulate values from logistic norm mix
simVals <- rlogitnorm(num,mu,sigma);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

