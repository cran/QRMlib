\name{rcopula.t}
\alias{rcopula.t}
\title{
  t Copula Simulation 
}
\description{
generates a random sample from the t copula 
}
\usage{
rcopula.t(n, df, Sigma=equicorr(d, rho), d=2, rho=0.7)
}
\arguments{
\item{n}{
number of observations 
}
\item{df}{
degrees of freedom 
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
This function is set up to allow quick simulation of t copulas 
with an equicorrelation structure. Simply enter a value for the 
dimension d and the correlation parameter rho. For more general 
correlation matrices specify Sigma. 
}
\seealso{
   \code{\link{rAC}},
   \code{\link{rcopula.gumbel}},
   \code{\link{rcopula.clayton}},
   \code{\link{rcopula.gauss}},
   \code{\link{rcopula.frank}}
}
\examples{
data <- rcopula.t(2000,df=4,d=6,rho=0.7); 
pairs(data); 
}
\keyword{distribution}

