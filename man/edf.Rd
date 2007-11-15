\name{edf}
\alias{edf}
\title{
  Empirical Distribution Function 
}
\description{
calculates the empirical distribution function at each element of a vector of observations 
}
\usage{
edf(v, adjust=FALSE)
}
\arguments{
\item{v}{
a vector 
}
\item{adjust}{
should the denominator be adjusted to be (n+1)? The default is FALSE 
}
}
\value{
vector of probabilities 
}
\examples{
data(smi);
data(ftse100);
TS1 <- window(ftse100, "1990-11-09", "2004-03-25");
TS1Augment <- alignDailySeries(TS1, method="before"); #gives 3490 observations
TS2Augment <- alignDailySeries(smi, method="before");
INDEXES.RAW <- merge(TS1Augment,TS2Augment);
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
plot(Udata);
}
\keyword{utilities}
