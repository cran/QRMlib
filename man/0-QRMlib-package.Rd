\name{QRMlib-package}
\alias{QRMlib-package}
\alias{QRMlib}
\docType{package}
\title{
This package provides R-language code to investigate concepts in a Quantitative Risk Management book
for those users without access to S-Plus.
}
\description{
This is a free R-language translation of the S-Plus library (QRMlib) designed to accompany the book 
\emph{Quantitative Risk Management: Concepts, Techniques and Tools}  by Alexander  J. McNeil, 
Rüdiger Frey and Paul Embrechts. It was built by Scott Ulman \email{scottulman@hotmail.com}.
A separate S-Plus version of the library can be downloaded from Alexander McNeil's URL.
}
\details{
\tabular{ll}{
Package: \tab QRMlib\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2007-04-08\cr
Depends: \tab R(>= 2.4.1),fCalendar,methods,fEcofin, mvtnorm, chron,its,Hmisc\cr
License: \tab GPL version 2 or newer\cr
URL: \tab http://www.math.etzh.ch/~mcneil/book/qrmlib.html\cr
Packaged: \tab April 8, 2007\cr
Built: \tab R 2.4.1; i386-pc-mingw32; 2007-04-08 12:00:00; windows\cr
}
The package provides an entire library of methods to investigate concepts associated
with QRM, including Market Risk, Credit Risk, and Operational Risk, as developed
in the textbook.  Additionally, it contains a set of chapter scripts which can be used
to build many of the graphs and tables in the text.  Under the library folder, look for
folders Chap2-Chap8 which contain the scripts.
}

\note{
The original S-Plus data files \bold{cac40, danish, DJ, dji, ftse100, FXGBP.RAW, hsi, 
nasdaq, nikkei, smi, sp500, xdax} are all S-Plus \sQuote{timeSeries} object files.  

Unfortunately, R-language has several different time-series classes, none of which coincides
with the S-Plus version.  The R-Metrics' class  \sQuote{timeSeries} (contained in library 
fCalendar) is the closest to an S-Plus timeSeries.  

Hence data files built in this R-language translation are all R-Metrics \sQuote{timeSeries}
types.  This means you must load the fCalendar library in order to use the R-data files.  
The setup program should require automatic loading of fCalendar when you open QRMlib. 

To automatically load the QRMlib package (and the fCalendar package), see 
\code{\link{profileLoadLibrary}}\cr
To automatically load the data files and save them in the current workspace, see
\code{\link{storeDataInWorkspace}}
}
\author{
S-Plus Original by Alexander McNeil; R-language port by Scott Ulman\cr

Maintainer: Scott Ulman <scottulman@hotmail.com> for R-language version
}
\references{
\emph{Quantitative Risk Management: Concepts, Techniques and Tools}  by Alexander  J. McNeil, Rüdiger Frey and Paul Embrechts\cr
Princeton Press, 2005
}
\seealso{
\code{\link{QRMBook-workspace}} and \code{\link{storeDataInWorkspace}} and 
\code{\link{profileLoadLibrary}}
}
\keyword{ package }

