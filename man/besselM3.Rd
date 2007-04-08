\name{besselM3}
\alias{besselM3}
\title{
  Modified Bessel Function of 3rd Kind 
}
\description{
calculates modified Bessel function of third kind 
}
\usage{
besselM3(lambda=9/2, x=2, logvalue=FALSE)
}
\arguments{
\item{lambda}{parameter of Bessel function}
\item{x}{2nd parameter of Bessel function}
\item{logvalue}{whether or not log value should be returned}
}
\value{
vector of values of Bessel function with same length as x 
}
\details{
see page 497 of QRM and references given there 
}
\seealso{
   \code{\link{rGIG}},
   \code{\link{dghyp}},
   \code{\link{dmghyp}}
}

\keyword{math}

