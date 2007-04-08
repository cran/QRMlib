\name{aggregateQuarterlySeries}
\alias{aggregateQuarterlySeries}
\title{aggregateQuarterlySeries() method }
\description{
This is one of several substitutes for the S-Plus language method\cr
\emph{aggregateSeries(timeseries, FUN={max, mean, colAvg,\cr
        colSums,...},by={weeks,months,quarters...},...)}.\cr 
The R-language \bold{aggregateQuarterlySeries()} function allows the user 
to calculate a less granular timeseries (monthly) from a daily time series 
by using a statistic like the max, mean, sum, colSums, etc. Note the R-methods 
do NOT contain a \bold{by=\dQuote{quarters}} parameter so the R-language user 
must select either the the \bold{aggregateWeeklySeries} method, the 
\bold{aggregateMonthlySeries()} method, or the 
\bold{aggregateQuarterlySeries()} method to get the desired result.
}
\details{ 
For example, the user might want to create a series of quarterly \bold{colSums} 
returns from a daily time series of returns.  Alternatively, (s)he 
might want the quarterly \bold{mean} series.  In either case, a less granular 
(quarterly) set of values is calculated from a daily timeSeries object.  Unfor-
tunately, the R-Metrics package has not yet implemented an R-version 
of the S-Plus aggregateSeries() method.

The \bold{aggregateWeeklySeries()}, \bold{aggregateMonthlySeries(),} and the
\bold{aggregateQuarterlySeries()} are interim functions developed to convert 
daily timeSeries to weekly, monthly, or quarterly timeSeries objects via a 
statistic like the max, mean, colAvg, or ColSums.


These functions exist in the functionsUtility.R file of the library.
}
\usage{
aggregateQuarterlySeries(timeseries, FUNC = colSums)
}

\arguments{
  \item{timeseries}{a (usually) daily timeSeries (R-Metrics type from 
   fCalendar) from which the user wants to extract a quarterly colSums 
   (or quarterly mean) timeSeries}
  \item{FUNC}{The name of a function to use in aggregating the data.  
             For example the max, mean, min, etc.  Default is 'colSums'.}
}

\value{
  A quarterly timeSeries object characterized by some statistic like mean, max, 
  min of the daily series over a quarter. The positions attribute 
  (dates <- rseries@positions ) of the new time series will be the LAST DAYS 
  OF THE RESPECTIVE QUARTERS for the timeSeries object.
}

\seealso{
   \code{\link{aggregateWeeklySeries}},
   \code{\link{aggregateMonthlySeries}}
}
\examples{
#load nasdaq data set:
data(nasdaq);
#Create daily return series:
nreturns <- -mk.returns(nasdaq) 
#convert to quarterly series using maximum value from each quarter:
quarterly.maxima <- aggregateQuarterlySeries(nreturns, FUNC=max);
\dontrun{
Ret.DJ <- mk.returns(DJ);
#Choose only 10 of the 30 stocks:
selection1 <- c("AXP","EK","BA","C","KO","MSFT","HWP",
                "INTC","JPM","DIS");
partialDJ30dailyTS <- Ret.DJ[,selection1];
partialDJ30daily <- cut(partialDJ30dailyTS, from="1992-12-31", 
                              to="2000-12-31");
partialDJ30quarterlyTS <- aggregateQuarterlySeries(partialDJ30daily, 
                           FUNC= colSums);
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{ts}

