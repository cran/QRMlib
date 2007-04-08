\name{ESnorm}
\alias{ESnorm}
\title{
  Expected Shortfall for Normal Distribution 
}
\description{
calculates expected shortfall for normal distribution 
}
\usage{
ESnorm(p, mean=0, sd=1)
}
\arguments{
\item{p}{
probability level 
}
\item{mean}{
mean 
}
\item{sd}{
standard deviation 
}
}
\value{
expected shortfall 
}
\details{
see page 45 of QRM 
}
\seealso{
   \code{\link{ESst}}
}

\examples{
ESnorm(c(0.95,0.99)) 
}
\keyword{distribution}

