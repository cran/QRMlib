\name{dji.df}
\alias{dji.df}
\docType{data}
\title{
  Dow Jones Index (dataframe Object) January 2, 1980-March 25, 2004.
  The .df indicates the dataframe object. 
}
\description{
The \code{dji.df} dataframe provides the daily closing value for the Dow Jones index 
from January 1980 to March 2004. QRMlib's R-version 1.4.2 and above supplies data 
in both timeSeries and data.frame versions.
}
\usage{
data(dji.df)
}
\format{
This dataframe object contains the prices for the index at dji.df[,2] and the corresponding dates 
at dji.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{dji}}
}
\keyword{datasets}
