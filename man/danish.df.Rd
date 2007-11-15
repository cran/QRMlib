\name{danish.df}
\alias{danish.df}
\title{
  Danish Data from January 1980 through December 1990 as data.frame object
}
\description{
The \code{danish.df} dataframe provides the daily closing value for the Danish fire losses 
in millions of kroner measured daily from January 1980 through December 1990. QRMlib's R-version 
1.4.2 and above supply data in both timeSeries and data.frame versions.
}
\usage{
data(danish.df)
}
\format{
This dataframe object contains the Fire losses (in million of kroner) at danish.df[,2] and the 
corresponding dates at danish.df[,1].  The dataframe can be converted to a timeSeries by calling 
the ConvertDFToTimeSeries() method in functionsUtility.R.
}
\seealso{
   \code{\link{danish}}
}
\keyword{datasets}


