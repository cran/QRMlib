\name{xdax.df}
\docType{data}
\alias{xdax.df}
\title{
  Xetra DAX German Index (timeSeries Object) January 3, 1994-March 25, 2004
}
\description{
The \code{xdax.df} dataframe provides the daily closing value for the German Xextra DAX index 
from January 1994 to March 2004. QRMlib's R-version 1.4.2 and above supplies data in both 
timeSeries and data.frame versions.
}
\usage{
data(xdax.df)
}
\format{
This dataframe object contains the prices for the index at xdax.df[,2] and the corresponding dates 
at xdax.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{xdax}}
}
\keyword{datasets}

