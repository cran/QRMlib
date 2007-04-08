\name{fit.norm}
\alias{fit.norm}
\title{
  Fit Multivariate Normal 
}
\description{
fits multivariate (or univariate) normal by maximum likelihood 
}
\usage{
fit.norm(data)
}
\arguments{
\item{data}{
matrix of data where each row is a vector observation 
}
}
\value{
list containing MLEs and value of likelihood at maximum 
}
\seealso{
   \code{\link{dmnorm}}
}
\examples{
data <- rmnorm(1000,rho=0.7,d=3); 
fit.norm(data); 
}
\keyword{methods}

