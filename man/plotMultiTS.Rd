\name{plotMultiTS}
\alias{plotMultiTS}
\title{ Plot Multiple Time Series}
\description{
Plots multiple timeSeries objects on the same graph.
Use this function if the plot.timeSeries() method from fCalendar returns the error
 \dQuote{Error in xy.coords(x, y, xlabel, ylabel, log) :} \cr 
   \dQuote{\sQuote{x} and sQuote{y} lengths differ}
}
\usage{
plotMultiTS(tS, colvec = 1:ncol(tS), type = "l", ltypvec = 1, lwdvec = 1, yrange, format, at, reference.grid, ...)
}
\arguments{
  \item{tS}{a timeSeries object with multiple data columns}
  \item{colvec}{By default all columns will be used.\cr  Use \emph{eqn{colvec = c(1,3)}} to use only columns 1 and 3}
  \item{type}{ the typical types for plot functions: defaults to \emph{lines} }
  \item{ltypvec}{the typical types for lines in plots; For example \eqn{1=solid} \cr and \eqn{2=dashed}. Will be expanded to a vector for columns selected. }
  \item{lwdvec}{ the typical line-width for lines in plots.}
  \item{yrange}{ the maximum and minimum values to appear on the y-axis of the plot }
  \item{format}{desired format for the date labels on x axis (defaults to year only) }
  \item{at}{the parameter for an \emph{axis.POSIXct} call.  Use only if you know the start and end dates 
and want a more granular set of labels on the x-axis than the default will provide. }
  \item{reference.grid}{ set to TRUE if you want vertical and horizontal grid lines added to plot }
  \item{\dots}{ Any additional parameters to pass to plot (such as \emph{par} parameters) }
}

\seealso{
   \code{\link{timeSeriesClass}}
}
\examples{
data(DJ);
Sdata <- window(DJ, from="1993-01-01", to="2000-12-31");
#select only 4 stocks from 30-stock index:
tsSelections <- c("GE","INTC","KO","JNJ");
Sdata <- Sdata[,tsSelections];
Zdata <- log(Sdata);
rm(Sdata);
#Plot all 4 columns on same graph:
plotMultiTS(Zdata, reference.grid=TRUE);
#plot only columns 2 and 3 on the graph:
plotMultiTS(Zdata, colvec= c(2,3),reference.grid=TRUE, format="\%Y-\%m");
}
\keyword{hplot }

