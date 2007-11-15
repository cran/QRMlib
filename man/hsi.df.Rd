\name{hsi.df}
\alias{hsi.df}
\docType{data}
\title{
  Hang Seng Stock Market Index (dataframe) January 1994 to March 2004
}
\description{
The \code{hsi.df} dataframe provides the daily closing value for the Hanh Seng Index 
from January 1994 to March 2004. QRMlib's R-version 1.4.2 and above 
supplies data in both timeSeries and data.frame versions.
}
\usage{
data(hsi.df)
}
\format{
This dataframe object contains the prices for the index at hsi.df[,2] and the corresponding dates 
at hsi.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{hsi}}
}
\keyword{datasets}

