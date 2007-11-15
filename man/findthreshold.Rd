\name{findthreshold}
\alias{findthreshold}
\title{
  Find a Threshold 
}
\description{
find threshold (or threshold vector) corresponding to given number of upper order statistics 
}
\usage{
findthreshold(data, ne)
}
\arguments{
\item{data}{
Data vector. See details section for extracting vector from other types.  
}
\item{ne}{
vector giving number of excesses above the threshold  
}
}
\value{
vector of suitable thresholds corresponding to each of the number of excesses given in the ne vector 
}
\details{
If using matrix as data, pass matname[,n] to pass nth column. \cr
If using a dataframe, pass dfname[[\dQuote{colname}]] or dfname[[n]] or dfname\$colname \cr
or dfname[ , \dQuote{colname}] or dfname[ ,n] where n is col number. \cr
If using a timeSeries, pass \dQuote{as.vector(tS@Data[,n])} to pass nth column of timeSeries data. \cr
When tied data values exist, a threshold is found so that at least 
the specified number of extremes lies above threshold. \cr
 
}
\seealso{
   \code{\link{fit.GPD}}
}

\examples{
#Load Danish data timeSeries file
data(danish);
targetVector <- as.vector(danish@Data);
# Find threshold giving (at least) fifty exceedances for Danish data 
findthreshold(targetVector,50);
}
\keyword{methods}

