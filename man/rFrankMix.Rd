\name{rFrankMix}
\alias{rFrankMix}
\title{
  Mixture Distribution Yielding Frank Copula 
}
\description{
generates random sample from discrete mixing distribution 
required for sampling from Frank's copula using Laplace 
transform method 
}
\usage{
rFrankMix(n, theta)
}
\arguments{
\item{n}{
size of sample 
}
\item{theta}{
value of parameter of Frank copula 
}
}
\value{
random sample of size n 
}
\details{
see page 224 of QRM. Algorithm generates V corresponding 
to Frank's copula. 
}
\seealso{
   \code{\link{rAC}}
}
\examples{
#Pass the parameter values n=20 and theta=0.5
result <- rFrankMix(20,0.5);
}
\keyword{distribution}

