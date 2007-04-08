\name{BiDensPlot}
\alias{BiDensPlot}
\title{
  Bivariate Density Plot 
}
\description{
makes perspective or contour plot of a bivariate density 
}
\usage{
BiDensPlot(func, xpts=c(-2, 2), ypts=c(-2, 2), npts=50, type="persp", ...)
}
\arguments{
\item{func}{
a function that evaluates on a n by 2 matrix to give n values of the bivariate density 
}
\item{xpts}{
limits of x range 
}
\item{ypts}{
limits of y range 
}
\item{npts}{
the number of subdivision points between x and y over the speicified range xpts to ypts
}
\item{type}{
"persp" or "contour" plot 
}
\item{\dots}{
further parameters of density function 
}
}
\section{Side Effects}{
produces a contour or perspective plot 
}
\seealso{
   \code{\link{dmnorm}},
   \code{\link{dmt}}
}
\examples{
BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7)) 
}
\keyword{hplot}
