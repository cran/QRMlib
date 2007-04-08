\name{timeSeriesClass}
\alias{timeSeriesClass}
\title{
  timeSeries Objects in R 
}
\description{
The R-language has developed multiple 'time-series' type objects 
through multiple contributors. The time-series objects in R-language 
closest to those in S-Plus appear to be those belonging to the timeSeries 
class described in fCalendar library from R-metrics.   Dates and times are
implemented as 'timeDate' objects wtihin 'timeSeries'.  The class contains 
functions for the generation and representation of 'timeSeries' objects and 
mathematical operations on the objects.\cr
Use timeSeries() as the constructor for the class.

Usage:\cr
timeSeries(data, charvec, units=NULL, format="ISO", zone="GMT",\cr 
     FinCenter=myFinCenter, \cr
     recordIDs=data.frame(),title=NULL, documentaion = NULL, ...)\cr  
}
\arguments{
\item{data}{
a vector or matrix or data frame containing numeric values 
}
\item{charvec}{
a character vector of dates and times 
}
\item{units}{
(optional) character string allowing overwrite of current column names  
}
\item{format}{
(optional) timeDate string format, defaulting to 'ISO' 
}
\item{zone}{
(optional)time zone where the data were recorded  
}
\item{FinCenter}{
location of financial center as Continent/City 
}
}
\value{
a timeSeries object 
}
\details{
IMPORTANT INFORMATION:
You can extract the DATE segment from a timeSeries object using\cr
1) the seriesPosition function:(e.g. use seriesPositions(sp500))\cr
2) the positions attribute (e.g. dates <- sp500@positions)
  
You can extract the NUMERIC segment from a timeSeries object using\cr 
1) the seriesData function (e.g. use seriesData(sp500))\cr
2) the Data attribute (e.g. returns <- sp500@Data) 
}
\seealso{
   \code{\link{TimeSeriesClassRMetrics}}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{ts}

