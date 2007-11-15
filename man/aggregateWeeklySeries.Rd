\name{aggregateWeeklySeries}
\alias{aggregateWeeklySeries}

\title{aggregateWeeklySeries() method }
\description{
This is one of several substitutes for the S-Plus language method\cr
\emph{aggregateSeries(timeseries, FUN={max, mean, colAvg,\cr
       colSums,...}, by={weeks,months,quarters...},...)}.\cr 
The R-language \bold{aggregateWeeklySeries()} function allows the user 
to calculate a less granular timeseries (weekly) from a daily time series 
by using a statistic like the max, mean, sum, etc. Note the R-methods do NOT 
contain a \bold{by=\dQuote{months}} parameter so the R-language user must select 
either the \bold{aggregateWeeklySeries} method, the \bold{aggregateMonthlySeries()} 
method, or the \bold{aggregateQuarterlySeries()} method to get the desired result.
}
\details{ 
For example, the user might want to create a series of weekly \bold{colSums} 
returns from a daily time series of returns.  Alternatively, (s)he 
might want the quarterly \bold{mean} series.  In either case, a less granular 
set of return values is calculated from a daily timeSeries object.  Unfortunately, 
the R-Metrics package has not yet implemented an R-version of the S-Plus 
aggregateSeries() method.

The \bold{aggregateWeeklySeries()}, \bold{aggregateMonthlySeries(),} and the
\bold{aggregateQuarterlySeries()} are interim functions developed to convert 
daily timeSeries to weekly, monthly, or quarterly timeSeries objects via a 
statistic like the max, mean, colAvg, or ColSums.


These functions exist in the functionsUtility.R file of the library.
}
\usage{
aggregateWeeklySeries(timeseries, FUNC = colSums)
}
\arguments{
  \item{timeseries}{a daily timeSeries (R-Metrics type from fCalendar) from 
which the user wants to extract a less granular weekly average timeSeries}
  \item{FUNC}{The name of a function to use in aggregating the data.  
             For example the colAvg, max, mean, min, etc.   The default is colSums}
}

\value{
  A weekly timeSeries object characterized by some statistic like colAvg, 
  of the daily series over a month. The positions attribute 
  (dates <- rseries@positions ) of the new time series will be the LAST DAYS 
  OF THE RESPECTIVE weeks for the timeSeries object.
}

\seealso{
   \code{\link{aggregateMonthlySeries}},
   \code{\link{aggregateQuarterlySeries}}
}
\examples{
#load nasdaq data set:
data(nasdaq);
data(DJ);
#Create daily return series:
nreturns <- -mk.returns(nasdaq) 
#convert to weekly series using colSums values (adding daily returns to get weekly)
weekly.nasdaq <- aggregateWeeklySeries(nreturns, FUNC=colSums);
Ret.DJ <- mk.returns(DJ);
#Choose only 10 of the 30 stocks:
selection1 <- c("AXP","EK","BA","C","KO","MSFT","HWP",
                "INTC","JPM","DIS");
partialDJ30dailyTS <- Ret.DJ[,selection1];
partialDJ30daily <- window(partialDJ30dailyTS, from="1993-01-01", 
                              to="2000-12-31");
partialDJ30weeklyTS <- aggregateWeeklySeries(partialDJ30daily, 
                       FUNC= colSums);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{ ts }
\keyword{ models}

