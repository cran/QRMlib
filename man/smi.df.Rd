\name{smi.df}
\docType{data}
\alias{smi.df}
\title{
  Swiss Market Index (dataframe Object) November 9, 1990 to March 25, 2004.
  The .df indicates the dataframe object.
}
\description{
The \code{smi.df} dataframe provides the daily closing value for the Swiss Market index 
from November 1990 to March 2004. QRMlib's R-version 1.4.2 and above supplies data in 
both timeSeries and data.frame versions.
}
\usage{
data(smi.df)
}
\format{
This dataframe object contains the prices for the index at smi.df[,2] and the corresponding dates 
at smi.df\$DATE. The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{smi}}
}
\keyword{datasets}

