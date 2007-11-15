\name{dcopula.t}
\alias{dcopula.t}
\title{
  t Copula Density 
}
\description{
evaluates density of t copula 
}
\usage{
dcopula.t(u, nu, P, logvalue=FALSE)
}
\arguments{
\item{u}{
matrix of dimension n times d, where d is the dimension of the copula 
and n is the number of vector values at which to evaluate density 
}
\item{nu}{
degrees of freedom of t copula 
}
\item{P}{
correlation matrix of t copula 
}
\item{logvalue}{
whether or not log density values should be returned (useful for ML) 
}
}
\value{
vector of density values of length n 
}
\details{
see pages 197 and 235 of QRM 
}
\seealso{
   \code{\link{dmt}},
   \code{\link{dcopula.clayton}},
   \code{\link{dcopula.gumbel}},
   \code{\link{dcopula.gauss}}
}
\examples{
ll <- c(0.01,0.99);
#create perspective plot for bivariate density:
BiDensPlot(func=dcopula.t,xpts=ll,ypts=ll,nu=4,P=equicorr(2,0.5));
}
\keyword{distribution}

