\name{sp500.df}
\docType{data}
\alias{sp500.df}
\title{
  Standard and Poors 500 Index (data.frame Object) January 2, 1990-March 25, 2004
}
\description{
The \code{sp500.df} data.frame provides the daily closing value for the S and P 500 Index 
from January 1980 to March 2004. QRMlib's R-version 1.4.2 and above supplies data in both 
timeSeries and data.frame versions.
}
\usage{
data(sp500.df)
}
\format{
This dataframe object contains the prices for the index at sp500.df[,2] and the corresponding dates 
at sp500.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.  
}
\seealso{
   \code{\link{sp500}}
}
\keyword{datasets}
