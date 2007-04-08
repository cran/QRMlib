\name{rtcopulamix}
\alias{rtcopulamix}
\title{
  Mixing Distribution on Unit Interval Yielding t Copula Model 
}
\description{
random generation for mixing distribution on unit interval 
yielding t copula model 
}
\usage{
rtcopulamix(n, pi, rho.asset, nu)
}
\arguments{
\item{n}{
sample size 
}
\item{pi}{
default probability 
}
\item{rho.asset}{
asset correlation parameter 
}
\item{nu}{
degree of freedom parameter 
}
}
\value{
random values on unit interval 
}
\details{
see page 361 in QRM; we consider exchangeable case of this model 
}
\seealso{
   \code{\link{rbeta}},
   \code{\link{rclaytonmix}},
   \code{\link{rlogitnorm}},
   \code{\link{rprobitnorm}}
}
\examples{
#set number, mean, and variance:
num <- 1000;
pi <- 0.9;
rho = 0.5;
df <- 4;
simVals <- rtcopulamix(num,pi,rho,df);
}
\keyword{distribution}

