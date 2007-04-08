\name{extremalPP}
\alias{extremalPP}
\title{Extremal Point Process }
\description{
  creates an extremal point process of class MPP
}
\usage{
extremalPP(data, threshold = NA, nextremes = NA)
}
\arguments{
  \item{data}{a timeSeries object or vector of numbers to be interpreted as a regular time series}
  \item{threshold}{threshold value (either this or "nextremes" must be given but not both)}
  \item{nextremes}{the number of upper extremes to be used (either this or "threshold" must 
 be given but not both)}
}
\details{
  see pages 298-301 of QRM 
}
\value{
  a list describing class MPP (marked point process) consisting of times and magnitudes 
   of threshold exceedances:
  \item{times}{vector of julian day counts (since 1/1/1960) for each exceedance}
  \item{marks}{vector of exceedances values (differences between value and threshold at each mark)}
  \item{startime}{the julian count one day prior to the first date in the entire timeSeries}
  \item{endtime}{value of last julian count in entire timeSeries}
  \item{threshold}{value of threshold above which exceedances are calculated}
}

\seealso{
   \code{\link{unmark}},
   \code{\link{fit.sePP}},
   \code{\link{fit.seMPP}}
}
\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500); 
tD <- timeDate("12/31/1995","\%m/\%d/\%Y"); 
window <- (seriesPositions(sp500.nreturns) > tD);
sp500.nreturns <- sp500.nreturns[window]; 
tmp <- extremalPP(sp500.nreturns,ne=100); 
tmp$marks[1:5];
tmp$threshold;
}
\keyword{models}

