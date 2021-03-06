\name{fit.gausscopula}
\alias{fit.gausscopula}
\title{
  Fit Gauss Copula 
}
\description{
fits Gauss copula to pseudo-copula data 
}
\usage{
fit.gausscopula(Udata)
}
\arguments{
\item{Udata}{
matrix of pseudo-copula data where rows are vector observations with all values in unit interval 
}
}
\value{
list containing parameter estimates and details of fit 
}
\details{
see pages 234-235 of QRM 
}
\seealso{
   \code{\link{fit.tcopula}},
   \code{\link{fit.Archcopula2d}}
}

\examples{
data(ftse100);
data(smi);
TS1 <- window(ftse100, "1990-11-09", "2004-03-25");
TS1Augment <- alignDailySeries(TS1, method="before");
TS2Augment <- alignDailySeries(smi, method="before");
INDEXES.RAW <- merge(TS1Augment,TS2Augment);
#Cleanup:
rm(TS1, TS1Augment, TS2Augment);
INDEXES <- mk.returns(INDEXES.RAW);
PARTIALINDEXES <- window(INDEXES, "1994-01-01", "2003-12-31");
#Now create a data matrix from the just-created timeSeries 
data <- seriesData(PARTIALINDEXES);
#Keep only the data items which are non-zero for both smi and ftse100
data <- data[data[,1]!=0 & data[,2] !=0,];
# Construct pseudo copula data. The 2nd parameter is MARGIN=2 
#when applying to columns and 1 applied to rows. Hence this says to
#apply the 'edf()' empirical distribtion function() to the columns
#of the data. 
Udata <- apply(data,2,edf,adjust=1);
mod.gauss <- fit.gausscopula(Udata); 
mod.gauss; 
}
\keyword{methods}

