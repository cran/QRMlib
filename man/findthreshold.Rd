\name{findthreshold}
\alias{findthreshold}
\title{
  Find a Threshold 
}
\description{
find threshold corresponding to given number of upper order statistics 
}
\usage{
findthreshold(data, ne)
}
\arguments{
\item{data}{
data vector 
}
\item{ne}{
vector giving number of extremes above the threshold  
}
}
\value{
vector of suitable thresholds 
}
\details{
when the data are tied a threshold is found so that at least 
the specified number of extremes lie above 
}
\seealso{
   \code{\link{fit.GPD}}
}

\examples{
data(danish);
# Find threshold giving (at least) fifty exceedances for Danish data 
findthreshold(danish,50);
}
\keyword{methods}

