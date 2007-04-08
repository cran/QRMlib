\name{mk.its.exceedances.tS}
\alias{mk.its.exceedances.tS}
\title{Extract Exceedances over a Threshold into Irregular Time Series for 
Plotting when Input is TimeSeries}
\description{
Builds an irregular time series object of the \sQuote{its} class 
to hold the exceedances above a specified threshold ordered 
according to date. The \sQuote{its} time series may then be plotted 
easily to show clustering effects (\dQuote{volatility clustering})
}
\details{
This function passes a timeSeries object containing both the 
Data and positions slots for associating dates with Data when 
building the irregular time series.   

See \emph{mk.its.exceedances.vector()} for a method which passes a 
separate vector of data and a matching-length position slot 
from a different timeSeries object. 
}
\usage{mk.its.exceedances.tS(timeseries, col = 1, thresholdValue)}
\arguments{
  \item{timeseries}{fCalendar timeSeries object containing data set
     for which you wish to plot exceedances over a threshold}
  \item{col}{the column number of the timeSeries data from which you 
             wish to extract the exceedances}
  \item{thresholdValue}{the minimum data value above which you want to 
   collect exceedances into the irregular time series (\sQuote{its} object); 
   calculate by calling QRM function findthreshold()}
}

\value{
An irregular time series (\sQuote{its} object requiring the \sQuote{its} package to 
be loaded) with exceedances ordered according to date.  Input \sQuote{its} 
into a plot() function to graph the exceedances.
}
\references{
see pages 117-8 of QRM.  Chapter 4 discusses the usefulness 
of graphing exceedances to display \dQuote{volatility clustering}.
}
\seealso{
   \code{\link{timeSeriesClass}},
   \code{\link{mk.its.exceedances.vector}},
   \code{\link{mk.returns}}
}

\examples{
data(sp500);
#Select data from the SP500 INDEX time series between 
#01/01/1998 and 12/31/2003 
index <- cut(sp500,"1998-01-01", "2003-12-31"); 
#Create a return time series using log differences. 
rseries <- mk.returns(index); 
kval <- 50; 
#set threshold so we have 50 observations above threshold 
#Find the threshold associated with the top 50 observations: 
upper <- findthreshold(rseries@Data ,kval) 
#Create an 'its' (irregular time series) associated with the 
#observations above the threshold 'upper'. 
itsSample <- mk.its.exceedances.tS(rseries,1,upper); 
#Plot the irregular time series of exceedances 
plot(itsSample,type="h", main="Clustering Evidence-S&P500 
       Index:1998-2003",ylab="Large Returns"); 
\dontrun{ 
#3rd stock in series of DJ 30 stocks from timeSeries DJ containing 
#all 30 stocks (not an index)
newSeries <- mk.returns(DJ); 
u3 <- findthreshold(newSeries@Data[,3],kval); 
#use 3rd column in Data slot 
itsSample2 <- mk.its.exceedances.tS(newSeries, 3, u3);  
#pass 3 for third column in Data slot 
plot(itsSample2,type="h", main="Clustering Evidence-DJ[,3]",
        ylab="Large Returns"); 
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{ts}
\keyword{methods}
