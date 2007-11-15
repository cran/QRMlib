\name{FXGBP.RAW.df}
\alias{FXGBP.RAW.df}
\docType{data}
\title{
  Sterling Exchange Rates as data.frame object January 1987 to March 2004.
  The .df indicates the dataframe object. 
}
\description{
The \code{FXGBP.RAW.df} dataframe provides daily exchange rates for major 
currencies (dollar, yen, euro, Swiss franc) against 
the British pound for the period January 1987 through March 2004 
}
\usage{
data(FXGBP.RAW.df)
}
\format{
This dataframe contains the following 5 columns: \cr
\tabular{ll}{
\emph{DATE} \tab the date for the corresponding rates \cr
\emph{GBP.USD} \tab exchange rate against U.S. dollar \cr
\emph{GBP.EUR} \tab exchange rate against Euro \cr
\emph{GBP.JYN} \tab exchange rate against Japanese yen   \cr
\emph{GBP.CHF} \tab exchange rate against Swiss Frank \cr
}
There are 4360 rows with daily rates from 1987 to 2004.
The dataframe can be converted to a timeSeries by calling the ConvertDFToTimeSeries()
method in functionsUtility.R.
}
\seealso{
   \code{\link{FXGBP.RAW}}
}
\keyword{datasets}
