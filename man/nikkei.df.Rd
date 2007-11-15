\name{nikkei.df}
\docType{data}
\alias{nikkei.df}
\title{
  Nikkei Stock Market Index (data.frame Object) January 4, 1994-March 25, 2004
}
\description{
The \code{nikkei.df} dataframe provides the daily closing value for the Nikkei index 
from January 1994 to March 2004 in its nikkei@Data slot. QRMlib's R-version 1.4.2 and above 
supplies data in both timeSeries and data.frame versions.
}
\usage{
data(nikkei.df)
}
\format{
This dataframe object contains the prices for the index at nikkei.df[,2] and the corresponding dates 
at nikkei.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{nikkei}}
}
\keyword{datasets}

