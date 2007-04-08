\name{fit.GEV}
\alias{fit.GEV}
\title{
  Fit Generalized Extreme Value Distribution  
}
\description{
fits generalized extreme value distribution (GEV) to block maxima data 
}
\usage{
fit.GEV(maxima)
}
\arguments{
\item{maxima}{
block maxima data 
}
}
\value{
list containing parameter estimates, standard errors and details of the fit 
}
\details{
see pages 271-272 of QRM 
}
\seealso{
   \code{\link{pGEV}},
   \code{\link{pGPD}},
   \code{\link{fit.GPD}}
}

\examples{
data(nasdaq);
nreturns <- -mk.returns(nasdaq); 
monthly.maxima <- aggregateMonthlySeries(nreturns,FUN=max); 
monthly.maxima <- seriesData(monthly.maxima) 
mod1 <- fit.GEV(monthly.maxima); 
}
\keyword{methods}

