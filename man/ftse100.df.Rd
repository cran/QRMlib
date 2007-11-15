\name{ftse100.df}
\alias{ftse100.df}
\docType{data}
\title{
  FTSE 100 Stock Market Index as dataframe object
}
\description{
The \code{ftse100} dataframe provides the daily closing value for the FTSE index 
from January 1980 to March 2004. QRMlib's R-version 1.4.2 and above supplies data 
in both timeSeries and data.frame versions.
}
\usage{
data(ftse100.df)
}
\format{
This dataframe object contains the prices for the index at ftse100.df[,2] and the corresponding dates 
at ftse100.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{ftse100}} 
}
\keyword{datasets}
