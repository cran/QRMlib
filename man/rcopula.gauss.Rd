\name{rcopula.gauss}
\alias{rcopula.gauss}
\title{
  Gauss Copula Simulation 
}
\description{
generates a random sample from the Gaussian copula 
}
\usage{
rcopula.gauss(n, Sigma=equicorr(d, rho), d=2, rho=0.7)
}
\arguments{
\item{n}{
number of observations 
}
\item{Sigma}{
correlation matrix 
}
\item{d}{
dimension of copula 
}
\item{rho}{
correlation parameter for specifying an equicorrelation structure 
}
}
\value{
a matrix with n rows and d columns 
}
\details{
This function is set up to allow quick simulation of Gauss copulas 
with an equicorrelation structure. Simply enter a value for the 
dimension d and the correlation parameter rho. For more general 
correlation matrices specify Sigma. 
}
\seealso{
   \code{\link{rAC}},
   \code{\link{rcopula.gumbel}},
   \code{\link{rcopula.clayton}},
   \code{\link{rcopula.frank}},
   \code{\link{rcopula.t}}
}
\examples{
data <- rcopula.gauss(2000,d=6,rho=0.7); 
pairs(data); 
}
\keyword{distribution}

