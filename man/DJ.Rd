\name{DJ}
\docType{data}
\alias{DJ}
\title{
  Dow Jones 30 Stock Prices (timeSeries object) January 1991 to December 2000
} 
\description{
The \code{DJ}  timeSeries data set provides the closing values of the Dow Jones 
30 Stocks from 1991-2000 in its DJ@Data slot. QRMlib's R-version 1.4.2 and above 
supplies data in both timeSeries and data.frame versions.
}
\usage{
data(DJ)
}
\format{
This timeSeries object contains the prices for the index at DJ@Data and the corresponding dates 
at DJ@positions. You may also access all prices of the first five stocks via DJ[,1:5].
}
\seealso{
   \code{\link{DJ.df}}
}
\keyword{datasets}


