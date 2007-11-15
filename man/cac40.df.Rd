\name{cac40.df}
\docType{data}
\alias{cac40.df}
\title{
  CAC 40 Stock Market Index (France) as dataframe object from anuary 1994 to March 25, 2004
}
\description{
The \code{cac40.df} dataframe data set provides the daily closing values of the French CAC 40 
stock index for the period 1994 to March 2004. QRMlib's R-version 1.4.2 and above 
supplies data in both timeSeries and data.frame versions.
}
\usage{
data(cac40)
}
\format{
This dataframe object contains the prices for the index at cac40.df[,2] and the corresponding dates 
at cac40.df[,1]. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R. 
}
\seealso{
   \code{\link{cac40}}
}
\keyword{datasets}
