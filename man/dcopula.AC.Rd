\name{dcopula.AC}
\alias{dcopula.AC}
\title{ Archimedean Copula Density }
\description{
  Evaluates the density associated with an Archimedean copula
}
\usage{
dcopula.AC(u, theta, name, logvalue = TRUE)
}
\arguments{
  \item{u}{matrix of dimension n times d, where d is the dimension of the copula and n is the number of 
vector values at which to evaluate density }
  \item{theta}{copula parameter }
  \item{name}{ copula name, e.g. "gumbel", "clayton" }
  \item{logvalue}{ whether or not log density values should be returned (useful for ML) }
}
\details{
 This is a generic function, designed so that further copulas, or expressions for densities of 
higher-dimensional copulas may be added. Clayton works in any dimension at present but Gumbel is only 
implemented for d=2. To extend one must calculate the dth derivative of the generator inverse and take 
logarithm of absolute value; this is the term called loggfunc. In addition, for other copulas, one 
needs the generator phi and the log of the negative value of its first derivative lnegphidash.
}
\value{
  vector of density values of length n
}
\seealso{ 
   \code{\link{dcopula.gauss}},
   \code{\link{fit.AC}}
}
\keyword{ distribution }

