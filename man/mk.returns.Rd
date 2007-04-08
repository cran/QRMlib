\name{mk.returns}
\alias{mk.returns}
\title{
  Make Financial Return Data 
}
\description{
makes financial return data from asset price data 
}
\usage{
mk.returns(tsdata, type="log")
}
\arguments{
\item{tsdata}{
a timeSeries object containing prices 
}
\item{type}{
whether "log" or "relative" returns should be constructed 
}
}
\value{
a timeSeries object containing returns 
}
\seealso{
   \code{\link{timeSeriesClass}},
   \code{\link{TimeSeriesClassRMetrics}}
}
\examples{
data(ftse100);
ftse100.r <- mk.returns(ftse100); 
}
\keyword{methods}

