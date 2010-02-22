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
Rudiger Frey and Paul Embrechts. It was built by Scott Ulman \email{scottulman@hotmail.com}.
A separate S-Plus version of the library can be downloaded from Alexander McNeil's URL.
}
\details{
\tabular{ll}{
Package: \tab QRMlib\cr
Type: \tab Package\cr
Version: \tab 1.4.5\cr
Date: \tab 2010-01-18\cr
Depends: \tab R(>= 2.7.0), methods, fSeries, mvtnorm, chron\cr
Suggests: its,Hmisc
License: \tab GPL version 2 or newer\cr
URL: \tab http://www.ma.hw.ac.uk/\eqn{\sim}{~}mcneil/book/index.html\cr
Packaged: \tab January 18, 2010\cr
Built: \tab R 2.10.1; i386-pc-mingw32; 2010-01-18 12:00:00; windows\cr
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
with the S-Plus version.  The R-Metrics' class \sQuote{timeSeries} (contained in package 
fSeries) is the closest to an S-Plus timeSeries. 

Unfortunately, RMetrics has made frequent significant changes to their timeSeries class.  
In fact, with R-2.6.0, RMetrics moved timeSeries from the fCalendar  
package to the fSeries package and made substantial changes.  This of course caused necessary
rewrites in QRMlib.

In R-2.10.0, RMetrics has moved the timeSeries object from fSeries package to a new timeSeries
package and changed the object.  Similarly, it moved the timeDate object from the fCalendar package
to a new timeDate package and changed the object.  

QRMlib continues to use the older timeDate object from fCalendar and the older timeSeries
object from fSeries.  It is \emph{very important} that you should NOT load the newer packages 
\sQuote{timeSeries} and \sQuote{timeDate} into a workspace where QRMlib is running.  
Loading those packages will \dQuote{mask} older functions with the same name and substitute 
their different behavior from the newer modules.  This will cause errors in some QRMlib functions.   

Although data files built in this R-language translation still use the older R-Metrics 
\sQuote{timeSeries} types from the fSeries package, I am now including another set of datasets 
built as dataframes which you may try to use in your analyses.\cr  
All data files are stored in the \emph{data} subfolder of QRMlib.  The timeSeries
file types end with the extension .rda while the dataframe types end with the extension .df.R.\cr 
See the section \emph{storeDataInWorkspace} for further details about the files.

To automatically load the QRMlib package, see 
\code{\link{profileLoadLibrary}}\cr
To automatically load the data files and save them in the current workspace, see
\code{\link{storeDataInWorkspace}}
}
\author{
S-Plus Original by Alexander McNeil; R-language port by Scott Ulman\cr

Maintainer: Scott Ulman <scottulman@hotmail.com> for R-language version
}
\references{
\emph{Quantitative Risk Management: Concepts, Techniques and Tools}  by Alexander  J. McNeil, Rudiger Frey and Paul Embrechts\cr
Princeton Press, 2005
}
\seealso{
\code{\link{QRMBook-workspace}},
\code{\link{storeDataInWorkspace}},  
\code{\link{profileLoadLibrary}}
}
\keyword{ package }

