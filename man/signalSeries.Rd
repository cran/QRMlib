\name{signalSeries}
\alias{signalSeries}
\title{signalSeries object}
\description{
  Structured after the S-Plus signalSeries object. It contains a 
  data slot of any type and a NUMERIC positions slot rather than 
  the date slot of a timeSeries. In other words, each data value 
  has a numeric value associated with its position in the overall list 
}
\usage{
signalSeries(data, positions., units, units.position, from = 1, by = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a component which is typically a dataframe}
  \item{positions.}{a numeric component describing the positions of the data values }
  \item{units}{ character vector describing the type of units used in the data structure }
  \item{units.position}{character vector describing the type of units used for the positions }
  \item{from}{ starting value of positions }
  \item{by}{ amount to skip between positions }
}
\details{
  If no arguments are supplied, the default (empty) signalSeries object is returned. Otherwise, 
a signalSeries object is created with the given positions and data, and units if they are 
supplied. As an alternative to supplying the positions directly, they can be supplied by giving 
from and by, in which case the positions are generated as a numeric sequence with the right length 
to match the data
}
\value{
  a signalSeries object with the given data and positions
}
\seealso{
   \code{\link{aggregateSignalSeries}}
}

\examples{
signalSeries(); #default object with no data or positions
#Create matrix of simulated values from multivariate-t distribution
m <- 90; n <- 3000;
dataSim <- rmt(m*n,df=3,rho=0.5,d=2);
dataSimSS <- signalSeries(dataSim); 
}
\keyword{ts}
\keyword{classes}

