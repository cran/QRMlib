\name{plot.MPP}
\alias{plot.MPP}
\title{ Plot Marked Point Process}
\description{
  creates a picture of a marked point process  
}

\usage{
\method{plot}{MPP}(x, \dots)}
\arguments{
\item{x}{ a point process object of class MPP (marked point process) }
\item{\dots}{
further parameters which may be passed to the plot function (see R help
about the plot function for further information)
}
}


\details{
  Creates an appropriate plot on graphical device. The input variable PP will be
internally separated into x and y values to pass to plot()  
}
\seealso{
   \code{\link{extremalPP}}
}
\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500);
window <- (seriesPositions(sp500.nreturns) > 
      timeDate("1995-12-31",format = "\%Y-\%m-\%d")) &
        (seriesPositions(sp500.nreturns) < 
            timeDate("2004-01-01",format = "\%Y-\%m-\%d"));
sp500.nreturns <- sp500.nreturns[window];
tmp <- extremalPP(sp500.nreturns,ne=100);
plot.MPP(tmp);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{hplot}

