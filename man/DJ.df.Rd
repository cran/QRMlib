\name{DJ.df}
\docType{data}
\alias{DJ.df}
\title{
  Dow Jones 30 Stock Prices (data.frame object) January 1991 to December 2000.
  The .df indicates the dataframe object. 
} 
\description{
The \code{DJ.df} dataframe provides the closing values of the Dow Jones 
30 Stocks from 1991-2000. QRMlib's R-version 1.4.2 and above supplies data in 
both timeSeries and data.frame versions.
}
\usage{
data(DJ.df)
}
\format{
This dataframe object contains the prices for all 30 stocks at DJ.df[,1:30] and the corresponding dates 
at DJ.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R
}
\seealso{
   \code{\link{DJ}}
}
\keyword{datasets}


