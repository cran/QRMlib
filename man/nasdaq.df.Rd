\name{nasdaq.df}
\docType{data}
\alias{nasdaq.df}
\title{
  NASDAQ Stock Market Index (data.frame object) January 3, 1994 to March 25, 2004
} 
\description{
The \code{nasdaq} timeSeries dataset provides the daily closing value for the daily 
closing values of the NASDAQ index from January 1994 to March 2004 in its nasdaq@Data slot. 
QRMlib's R-version 1.4.2 and above  supplies data in both timeSeries and data.frame versions.
}
\usage{
data(nasdaq)
}
\format{
This dataframe object contains the prices for the index at nasdaq.df[,2] and the corresponding dates 
at nasdaq.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{nasdaq}}
}
\keyword{datasets}

