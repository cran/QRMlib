\name{plot.sePP}
\alias{plot.sePP}
\title{ Plot Self-Exciting Point Process}
\description{
 plots a fitted self-exciting point process model, either unmarked or marked 
}
\usage{
\method{plot}{sePP}(x,\dots)}
\arguments{
\item{x}{a fitted self-exciting point process model created by either 
fit.sePP or fit.seMPP. \sQuote{x} is the generic value passed to all S3
plot functions. }
\item{\dots}{
further parameters which may be passed to the plot function (see R help
about the plot function for further information)
}
}
\details{
Creates an appropriate plot on graphical device. The input variable will be
internally separated into x and y values to pass to plot()  
}
\seealso{
   \code{\link{fit.sePP}},
   \code{\link{fit.seMPP}}
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
mod2a <- fit.sePP(tmp,mark.influence=FALSE,std.errs=TRUE); 
plot.sePP(mod2a);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{hplot}
