\name{dcopula.clayton}
\alias{dcopula.clayton}
\title{
  Bivariate Clayton Copula Density 
}
\description{
evaluates density of bivariate Clayton copula 
}
\usage{
dcopula.clayton(u, theta, logvalue=FALSE)
}
\arguments{
\item{u}{
matrix of dimension n times 2, where 2 is the dimension of the copula and 
n is the number of vector values at which to evaluate density 
}
\item{theta}{
parameter of Clayton copula 
}
\item{logvalue}{
whether or not log density values should be returned (useful for ML) 
}
}
\value{
vector of density values of length n 
}
\details{
see page 192 of QRM for Clayton copula 
}
\seealso{
   \code{\link{fit.Archcopula2d}},
   \code{\link{dcopula.gauss}},
   \code{\link{dcopula.t}},
   \code{\link{dcopula.gumbel}}
}
\examples{
\dontrun{
#define a function:
normal.metaclayton <- function(x,theta);
{
  exp(dcopula.clayton(apply(x,2,pnorm),theta,logvalue=TRUE) + 
  apply(log(apply(x,2,dnorm)),1,sum));
}
#use function to create perspective plot for bivariate density: 
BiDensPlot(normal.metaclayton,xpts=ll,ypts=ll,npts=80,theta=2.2);
}
}
\keyword{distribution}
