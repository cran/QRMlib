\name{spdata.raw}
\alias{spdata.raw}
\docType{data}
\title{
  Standard and Poors Default Data 
}
\description{
The \code{spdata.raw} data frame has 20 rows and 11 columns. 
It contains default data for A, BBB, BB, B and C-rated companies for the years 1981 to 2000 
}
\usage{
data(spdata.raw)
}
\format{
This data frame contains the following 11 columns: \cr
\tabular{ll}{
\emph{year} \tab year of default\cr
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
   \code{\link{spdata}},
   \code{\link{momest}}
}
\examples{
data(spdata.raw);
attach(spdata.raw); 
BdefaultRate <- Bdefaults/Bobligors;
BBdefaultRate <- BBdefaults/BBobligors;
BBBdefaultRate <- BBBdefaults/BBBobligors; 
AdefaultRate <- Adefaults/Aobligors; 
CCCdefaultRate <- CCCdefaults/CCCobligors; 
\dontrun{
plot(year,CCCdefaultRate,xlab="Year",ylab="Rate",type="l"); 
lines(year,BdefaultRate,col=2); 
lines(year,BBdefaultRate,col=3); 
lines(year,BBBdefaultRate,col=4); 
lines(year,AdefaultRate,col=5); 
momest(Bdefaults,Bobligors); 
pi.B <- momest(Bdefaults, Bobligors)[1]; 
pi2.B <- momest(Bdefaults, Bobligors)[2];
rhoY.B <- (pi2.B-pi.B^2)/(pi.B-pi.B^2); 
}
}
\keyword{datasets}

