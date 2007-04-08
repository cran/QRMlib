\name{plot.timeSeriesIts}
\alias{plot.timeSeriesIts}
\title{ Plot single or multiple timeSeries objects on same graph }
\description{
Plots individual or multiple timeSeries objects with time on the x-axis 
and corresponding sereies values on the y-axis. Builds a reference grid 
by default. This function converts the timeSeries input object to an 
\sQuote{its}(irregular time series object) and hence requires the 
\sQuote{its} library to be loaded.\cr
Use plot.timeSeriesIts()whenever the timeSeries x has multiple columns 
of data.
}
\usage{
plot.timeSeriesIts(x, reference.grid = TRUE, lty = 1, ...)
}
\arguments{
  \item{x}{ a timeSeries object (may have multiple series) }
  \item{reference.grid}{ TRUE if you desire the graph to have a grid }
  \item{lty}{ line type for graphed series }
  \item{\dots}{ any other plot attributes which should be passed }
}
\details{
 plot.timeSeries() is not working properly in fCalendar 240.10068 for
 multiple timeSeries on the same graph.  Hence we've reinstituted an older
 version and named it plot.timeSeriesIts() to indicate we must have the
 \sQuote{its} library loaded.  We may deprecate this function eventually.
}
\author{documentation by Scott Ulman for R-language distribution}
\examples{
#Plot two time series (ftse100 and smi indices) on same graph
data(ftse100);
data(smi);
TS1 <- cut(ftse100, "1990-11-08", "2004-03-25");
TS1Augment <- alignDailySeries(TS1, method="before");
TS2Augment <- alignDailySeries(smi, method="before");
INDEXES.RAW <- merge(TS1Augment,TS2Augment);
#Cleanup:
rm(TS1, TS1Augment, TS2Augment);
plot.timeSeriesIts(INDEXES.RAW);
}

\keyword{hplot }

