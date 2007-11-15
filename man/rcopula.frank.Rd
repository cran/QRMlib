\name{rcopula.frank}
\alias{rcopula.frank}
\title{
  Frank Copula Simulation 
}
\description{
generates a random sample from the Frank copula 
}
\usage{
rcopula.frank(n, theta, d)
}
\arguments{
\item{n}{
sample size 
}
\item{theta}{
parameter value 
}
\item{d}{
dimension of copula 
}
}
\value{
matrix with n rows and d columns where rows are realizations 
}
\details{
see pages 222-224 in QRM 
}
\seealso{
   \code{\link{rAC}},
   \code{\link{rcopula.gumbel}},
   \code{\link{rcopula.clayton}},
   \code{\link{rcopula.gauss}},
   \code{\link{rcopula.t}}
}
\examples{
#generate data for a 2-dimensional frank copula with theta=4 
simdata <- rcopula.frank(1000,2,4); 
pairs(simdata); 
simdata <- rAC("frank",1000,2,4)
}
\keyword{distribution}

