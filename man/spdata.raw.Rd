\name{spdata.raw}
\alias{spdata.raw}
\docType{data}
\title{
  Standard and Poors Default Data (timeSeries object)
}
\description{
The \code{spdata.raw} timeSeries has 20 rows (dates) and 10 columns (obligors). 
It contains default data for A, BBB, BB, B and C-rated companies for the years 1981 to 2000. 
}
\usage{
data(spdata.raw)
}
\format{
This timeSeries contains the following 10 columns: \cr
\tabular{ll}{
\emph{Aobligors} \tab number of A-rated companies \cr
\emph{Adefaults} \tab number of A-rated companies defaulting in year \cr
\emph{BBBobligors} \tab number of BBB-rated companies  \cr
\emph{BBBdefaults} \tab number of BBB-rated companies that default in year \cr
\emph{BBobligors} \tab number of BB-rated companies  \cr
\emph{BBdefaults} \tab number of BB-rated companies that default in year \cr
\emph{Bobligors} \tab number of B-rated companies  \cr
\emph{Bdefaults} \tab number of B-rated companies that default in year \cr
\emph{CCCobligors} \tab number of CCC-rated companies  \cr
\emph{CCCdefaults} \tab number of CCC-rated companies that default in year \cr
}
There are 20 rows with values for the years from 1981 to 2000
}
\source{
Standard & Poors Credit Monitor 
}
\seealso{
   \code{\link{spdata.raw.df}},
   \code{\link{spdata}},
   \code{\link{momest}}
}
\examples{
data(spdata.raw);
attach(spdata.raw); 
BdefaultRate <- Bdefaults/Bobligors;
#Get an array of BB default rates for the 20 years 1981-2000:
BBdefaultRate <- BBdefaults/BBobligors;
BBBdefaultRate <- BBBdefaults/BBBobligors; 
AdefaultRate <- Adefaults/Aobligors; 
CCCdefaultRate <- CCCdefaults/CCCobligors; 
}
\keyword{datasets}

