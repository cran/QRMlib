\name{rcopula.clayton}
\alias{rcopula.clayton}
\title{
  Clayton Copula Simulation 
}
\description{
generates a random sample from the Clayton copula 
}
\usage{
rcopula.clayton(n, theta, d)
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
   \code{\link{rcopula.gauss}},
   \code{\link{rcopula.t}},
   \code{\link{rcopula.frank}}
}
\examples{
data <- rcopula.clayton(1000,2,4); 
pairs(data); 
}
\keyword{distribution}

