\name{equicorr}
\alias{equicorr}
\title{
  Equicorrelation Matrix 
}
\description{
constructs an equicorrelation matrix 
}
\usage{
equicorr(d, rho)
}
\arguments{
\item{d}{
dimension of matrix 
}
\item{rho}{
value of correlation 
}
}
\value{
an equicorrelation matrix 
}
\seealso{
   \code{\link{rmnorm}},
   \code{\link{rmt}}
}

\examples{
equicorr(7,0.5);
# Bivariate Visualization
ll <- c(0.01,0.99)
BiDensPlot(func=dcopula.gauss,xpts=ll,ypts=ll,P=equicorr(2,0.5));
BiDensPlot(func=dcopula.t,xpts=ll,ypts=ll,nu=4,P=equicorr(2,0.5)); 
}

\keyword{array}

