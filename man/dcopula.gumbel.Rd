\name{dcopula.gumbel}
\alias{dcopula.gumbel}
\title{
  Bivariate Gumbel Copula Density}
\description{
evaluates density of bivariate Gumbel copula}
\usage{
dcopula.gumbel(u, theta, logvalue=FALSE)
}
\arguments{
\item{u}{
matrix of dimension n times 2, where 2 is the dimension of the copula and n 
is the number of vector values at which to evaluate density 
}
\item{theta}{parameter of Gumbel copula}
\item{logvalue}{whether or not log density values should be returned 
(useful for ML)}
}
\value{
vector of density values of length n }
\details{
see page 192 of QRM for Gumbel copula }
\seealso{
   \code{\link{fit.Archcopula2d}},
   \code{\link{dcopula.clayton}},
   \code{\link{dcopula.t}},
   \code{\link{dcopula.gauss}}
}
\examples{
\dontrun{
normal.metagumbel <- function(x,theta)
{
  exp(dcopula.gumbel(apply(x,2,pnorm),theta,logvalue=TRUE) + 
  apply(log(apply(x,2,dnorm)),1,sum));
}
#use function to create perspective plot for bivariate density: 
BiDensPlot(normal.metagumbel,xpts=ll,ypts=ll,npts=80,theta=2);
}
}
\keyword{distribution}

