\name{dcopula.gauss}
\alias{dcopula.gauss}
\title{
  Gauss Copula Density 
}
\description{
evaluates density of Gauss copula 
}
\usage{
dcopula.gauss(u, P, logvalue=FALSE)
}
\arguments{
\item{u}{
matrix of dimension n times d, where d is the dimension of the copula 
and n is the number of vector values at which to evaluate density 
}
\item{P}{
correlation matrix of Gauss copula 
}
\item{logvalue}{
whether or not log density values should be returned (useful for ML) 
}
}
\value{
vector of density values of length n 
}
\details{
see pages 197 and 234 in QRM 
}
\seealso{
   \code{\link{dmnorm}},
   \code{\link{dcopula.clayton}},
   \code{\link{dcopula.t}},
   \code{\link{dcopula.gumbel}}
}
\examples{
ll <- c(0.01,0.99);
#create perspective plot for bivariate density:
BiDensPlot(func=dcopula.gauss,xpts=ll,ypts=ll,P=equicorr(2,0.5));
}
\keyword{distribution}

