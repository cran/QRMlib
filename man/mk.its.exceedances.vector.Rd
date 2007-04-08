\name{mk.its.exceedances.vector}
\alias{mk.its.exceedances.vector}
\title{Extract Exceedances over Threshold into Irregular Time Series for 
Plotting when Input is Vector}
\description{
Builds an irregular time series object of the \sQuote{its} class to hold 
the exceedances above a specified threshold ordered according to date. 
The \sQuote{its} time series may then be plotted easily to show clustering 
effects (\dQuote{volatility clustering}). 
}
\details{
This function passes a vector (not a matrix or timeSeries) containing the data.  
It passes a separate positions slot from a timeSeries for associating dates with 
Data when building the irregular time series.  

See \emph{mk.its.exceedances.tS()} for a method which passes a 
timeSeries object containing both data and position.  
}
\usage{mk.its.exceedances.vector(datavector, paralleltimeseriesPos, thresholdValue)}
\arguments{
  \item{datavector}{a data vector for which you wish to plot exceedances over a threshold}
  \item{paralleltimeseriesPos}{ a vector of length equal to the datavector containing 
    the positions (time-date) attribute from an R-Metrics type timeSeries object} 
  \item{thresholdValue}{the minimum data value above which you want to collect 
     exceedances into the irregular time series (\sQuote{its} object); 
     calculate by calling QRM function findthreshold()}
}

\value{
An irregular time series \sQuote{its} object requiring the \sQuote{its} package to 
be loaded with exceedances ordered according to date.  Input \sQuote{its} 
into plot function to graph the exceedances
}
\references{
See pages 117-8 of QRM.  Chapter 4 discusses the usefulness 
of graphing exceedances to display \dQuote{volatility clustering}.
}
\seealso{
   \code{\link{timeSeriesClass}},
   \code{\link{mk.its.exceedances.tS}},
   \code{\link{mk.returns}}
}
\examples{
\dontrun{
#In this example, a simulated data series tdata was 
#generated using parameters derived from fitting a timeSeries 
#for the SP500 index called rseries over a specified period. 
#tdata was randomly generated from the t-distribution using 
#parameters estimated from the actual SP500 time series of data 
#returns (log differences). Thus tdata is a VECTOR (one column) 
#rather than a timeSeries or matrix.  If you have a matrix, pass 
#only the single column you want to analyze. 
#set threshold so we have 50 observations above threshold from tdata.
kval <- 50;
upper <- findthreshold(tdata,kval); 
#Create an its (irregular time series) associated with observations 
#above threshold upper. We pass the simulated data vector tdata, 
#the dates (positions slot) from a separate timeseries over the same 
#period for which we are considering the simulated data, 
#and the threshold value for our data vector. 
itsSample <- mk.its.exceedances.vector(tdata,rseries@positions,upper); 
#Plot the irregular time series of exceedances 
plot(itsSample,type="h", main="Clustering Evidence-Simulated t data",
   ylab="Large Returns"); 
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}
\keyword{ts}

