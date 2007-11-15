\name{plot.PP}
\alias{plot.PP}
\title{ Plot Point Process}
\description{
  creates a picture of an unmarked point process.
}
\usage{
\method{plot}{PP}(x, \dots)}
\arguments{
\item{x}{a point process object of class PP which must be unmarked}
\item{\dots}{
further parameters which may be passed to the plot function (see R help
about the plot function for further information)
}
}
\details{
Creates an appropriate plot on graphical device. The input variable x will be
internally separated into starttime and endtime values to pass to plot.stepfun()    
}
\seealso{
   \code{\link{extremalPP}},
   \code{\link{unmark}}
}

\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500);
window <- (seriesPositions(sp500.nreturns) > 
      timeDate("1995-12-31",format = "\%Y-\%m-\%d")) &
        (seriesPositions(sp500.nreturns) < 
            timeDate("2004-01-01",format = "\%Y-\%m-\%d"));
sp500.nreturns <- sp500.nreturns[window];
#The following functions are contained in functionsHawkes.R. 
#Plot the 100 largest exceedances. Create an MPP (marked point process) class 
tmp <- extremalPP(sp500.nreturns,ne=100);
#Be sure to graph with plot.PP instead of plot.MPP:
tmp2 <- unmark(tmp);
plot.PP(tmp2);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{ hplot }
