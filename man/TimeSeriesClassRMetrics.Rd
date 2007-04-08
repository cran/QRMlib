\name{TimeSeriesClassRMetrics}
\alias{TimeSeriesClassRMetrics}
\title{
  timeSeries Class and Methods 
}
\description{
A collection and description of functions and methods dealing 
with regular and irregular 'timeSeries' objects. Dates and times 
are implemented as 'timeDate' objects. Included are functions and 
methods for the generation and representation of 'timeSeries' 
objects, and for mathematical operations. 

The functions and methods for the Generation of 'timeSeries' 
Objects are:
\tabular{ll}{
\sQuote{timeSeries()} \tab Creates a 'timeSeries' object from scratch,\cr  
\sQuote{read.timeSeries()} \tab Reads a \sQuote{timeSeries} from a spreadsheet file,\cr  
\sQuote{as.timeSeries()} \tab S3: Creates \sQuote{time Series} from a \sQuote{matrix},\cr  
\sQuote{is.timeSeries()} \tab S3: Tests if an object is of class a \sQuote{timeSeries},\cr  
\sQuote{print.timeSeries()} \tab S3: Print method for a \sQuote{timeSeries} object,\cr  
\sQuote{plot.timeSeries()} \tab S3: Plot method for a \sQuote{timeSeries} object,\cr  
\sQuote{lines.timeSeries()} \tab S3: Lines method for a \sQuote{timeSeries} object,\cr  
\sQuote{Ops.timeSeries()} \tab S3: Arith method for a \sQuote{timeSeries} object,\cr  
\sQuote{[.timeSeries()} \tab S3: \dQuote{[} method for a \sQuote{timeSeries} object,\cr  
\sQuote{head.timeSeries()} \tab S3: returns the head of a \sQuote{timeSeries} object,\cr  
\sQuote{tail.timeSeries()} \tab S3: returns the tail of a \sQuote{timeSeries} object.
}

The functions and methods for the Representation of \sQuote{timeSeries} 
Objects are:
\tabular{ll}{
\sQuote{seriesData()} \tab  Extracts data slot from a 'timeSeries',\cr  
\sQuote{seriesPositions()} \tab Extracts positions slot from a 'timeSeries',\cr  
\sQuote{start.timeSeries()} \tab S3: Extracts start date of a 'timeSeries',\cr  
\sQuote{end.timeSeries()} \tab S3: Extracts end date of a 'timeSeries',\cr  
\sQuote{as.vector.timeSeries()} \tab S3: Converts a 'timeSeries' to a vector,\cr  
\sQuote{as.matrix.timeSeries()} \tab S3: Converts a 'timeSeries' to a matrix,\cr  
\sQuote{as.data.frame.timeSeries()} \tab S3: Converts a 'timeSeries' to a data.frame.
}

The functions and methods for Math Operations of 'timeSeries' 
Objects are: 
\tabular{ll}{
\sQuote{applySeries()} \tab Applies a function to margins of a \sQuote{timeSeries},\cr  
\sQuote{alignDailySeries()} \tab Aligns a daily \sQuote{timeSeries} to new positions,\cr  
\sQuote{cut()} \tab Cuts out a piece from a \sQuote{timeSeries} object,\cr  
\sQuote{mergeSeries()} \tab Merges a \sQuote{timeSeries} object with a \sQuote{matrix},\cr  
\sQuote{ohlcDailyPlot()} \tab Plots open high low close bar chart,\cr  
\sQuote{revSeries()} \tab Reverts the order of \sQuote{timeSeries} object,\cr  
\sQuote{diffSeries()} \tab Takes differences from a \sQuote{timeSeries} object,\cr  
\sQuote{lagSeries()} \tab Lags a \sQuote{timeSeries} object,\cr  
\sQuote{outlierSeries()} \tab Removes outliers from a \sQuote{timeSeries} object,\cr  
\sQuote{returnSeries()} \tab Computes returns from a \sQuote{timeSeries} object,\cr  
\sQuote{logSeries()} \tab Returns logarithms of a \sQuote{timeSeries} object,\cr  
\sQuote{absSeries()} \tab Returns absolute values of a \sQuote{timeSeries} object. 
}


Functions calls include:\cr
timeSeries(data, charvec, units = NULL, format = "ISO", zone = "GMT", 
    FinCenter = myFinCenter, recordIDs = data.frame(), title = NULL, 
    documentation = NULL, ...)\cr 

read.timeSeries(file, zone = "GMT", FinCenter = "", title = "", 
    documentation = "", sep = ";")\cr
    
as.timeSeries(x, dimnames = TRUE, format = "")\cr
is.timeSeries(object)\cr 
The following are S3 method for class 'timeSeries':\cr
print(x, ...)\cr
plot(x, reference.grid = TRUE, lty = 1, ...)\cr 
lines(x, ...)\cr
Ops(e1, e2)\cr
x[i = min(1, nrow(x@Data]):nrow(x@Data),\cr 
  j = min(1, ncol(x@Data)):ncol(x@Data))\cr
head(x, ...)\cr
tail(x, ...)\cr

seriesData(object)\cr
seriesPositions(object)\cr
The following are S3 method for class 'timeSeries':\cr
start(x, ...)\cr
end(x, ...)\cr
as.vector(x, mode = "any")\cr
as.matrix(x)\cr
as.data.frame(x, row.names = NULL, optional = NULL)\cr

applySeries(x, from = NULL, to = NULL, by=c("monthly","quarterly"), FUN = colAvgs, units = NULL, ...);\cr 
cut(x, from, to);\cr
diffSeries(x, lag = 1, diff = 1, trim = FALSE, pad = NA);\cr
lagSeries(x, k = 1, trim = FALSE, units = NULL);\cr
outlierSeries(x, sd = 10, complement = TRUE);\cr
mergeSeries(x, y, units = NULL);\cr
returnSeries(x, type = c("continuous", "discrete"), percentage = FALSE, 
    trim = TRUE, digits = 4, units = NULL);\cr
revSeries(x);\cr
logSeries(x);\cr
absSeries(x);\cr
alignDailySeries(x, method = c("before", "after", "interp", "fillNA"), 
    include.weekends = FALSE, units = NULL);\cr
ohlcDailyPlot(x, volume = TRUE, colOrder = c(1:5), units = 1e6, 
    xlab = c("Date", "Date"), ylab = c("Price", "Volume"), 
    main = c("O-H-L-C", "Volume"), grid.nx = 7, grid.lty = "solid", ...);\cr 
}


\details{
Generation of Time Series Objects:\cr 
We have defined a timeSeries class which is in many aspects 
similar to the S-Plus class with the same name, but has also 
some important differences. The class has seven Slots, the 
\sQuote{Data} slot which holds the time series data in matrix form, 
the \sQuote{position} slot which holds the time/date as a character 
vector, the \sQuote{format} and \sQuote{FinCenter} slots which are the same 
as for the \sQuote{timeDate} object, the \sQuote{units} slot which holds 
the column names of the data matrix, and a \sQuote{title} and a 
\sQuote{documentation} slot which hold descriptive character strings. 
Date and time is managed in the same way as for timeDate objects. 
 
}
\seealso{
   \code{\link{timeSeriesClass}}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{ts}
\keyword{classes}
