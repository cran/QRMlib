\name{hessb}
\alias{hessb}
\title{
  Approximate Hessian Matrix 
}
\description{
calculates a numerical approximation of Hessian matrix 
}
\usage{
hessb(f, x, ep=0.0001, ...)
}
\arguments{
\item{f}{
function 
}
\item{x}{
value of function at which to approximate Hessian 
}
\item{ep}{
precision for numerical differencing 
}
\item{\dots}{
other arguments of function f 
}
}
\value{
matrix of approximate second derivatives 
}
\examples{
\dontrun{
#within fit.NH we approximate 2nd derivatives to calc standard errors
if(se) 
{
  hessmatrix <- hessb(negloglik,par.ests)
  vcmatrix <- solve(hessmatrix)
  par.ses <- sqrt(diag(vcmatrix))
  names(par.ses) <- names(par.ests)
  dimnames(vcmatrix) <- list(names(par.ests), names(par.ests))
}
else 
{
  par.ses <- NA
  vcmatrix <- NA
}
}
}
\keyword{utilities}

