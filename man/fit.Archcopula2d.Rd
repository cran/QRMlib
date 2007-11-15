\name{fit.Archcopula2d}
\alias{fit.Archcopula2d}
\title{
  Fit 2D Archimedean Copula 
}
\description{
Fits two-dimensional Archimedean copula by maximum likelihood.
This function has been deprecated.  Use fit.AC() instead. 
}
\usage{
fit.Archcopula2d(Udata, name)
}
\arguments{
\item{Udata}{
Matrix of copula data with two columns taking values in unit interval (hence Udata). 
}
\item{name}{
name of Archimedean copula: "clayton", "gumbel" 
}
}
\value{
list containing parameter estimate, standard error, value of log-likelihood 
at maximum and convergence flag 
}
\details{
see pages 234-236 of QRM 
}
\seealso{
   \code{\link{fit.gausscopula}},
   \code{\link{fit.tcopula}}
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
Udata <- apply(data,2,edf,adjust=TRUE);
#Fit 2-dimensional Archimedian copula: choices are gumbel or clayton
#using pseudo data generated via edf() from observed data:
mod.gumbel <- fit.Archcopula2d(Udata,"gumbel");
mod.clayton <- fit.Archcopula2d(Udata,"clayton");
}
\keyword{methods}

