\name{rcopula.gumbel}
\alias{rcopula.gumbel}
\title{
  Gumbel Copula Simulation 
}
\description{
generates a random sample from the Gumbel copula 
}
\usage{
rcopula.gumbel(n, theta, d)
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
   \code{\link{rcopula.frank}}
   \code{\link{rcopula.clayton}},
   \code{\link{rcopula.gauss}},
   \code{\link{rcopula.t}}
}
\examples{
data <- rcopula.gumbel(1000,3,4); 
pairs(data); 
}
\keyword{distribution}

