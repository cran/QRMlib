\name{plotTail}
\alias{plotTail}
\title{
  Tail Plot of GPD Model 
}
\description{
plots the tail estimate corresponding to a GPD model 
of excesses over a high threshold 
}
\usage{
plotTail(object, extend=2, fineness=1000, ...)
}
\arguments{
\item{object}{
result of fitting GPD to excesses over a high threshold 
}
\item{extend}{
how far plot should extend expressed as multiple of largest data value 
}
\item{fineness}{
number of points at which to evaluate the tail estimate 
}
\item{\dots}{
additional arguments for plot function 
}
}
\section{Side Effects}{
a plot of the tail estimate is produced on a graphical device 
}
\details{
see pages 282-284 in QRM 
}
\seealso{
   \code{\link{fit.GPD}},
   \code{\link{MEplot}}
}
\examples{
data(danish);
mod <- fit.GPD(danish,10); 
mod$par.ests; 
plotTail(mod); 
}
\keyword{hplot}

