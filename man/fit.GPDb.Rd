\name{fit.GPDb}
\alias{fit.GPDb}
\title{
  Fit Generalized Pareto Model B 
}
\description{
fits a generalized Pareto distribution to threshold exceedances using nlminb()
rather than nlmin()
}
\usage{
fit.GPDb(data, threshold=NA, nextremes=NA, method="ml", information="observed")
}
\arguments{
\item{data}{
data vector or times series 
}
\item{threshold}{
a threshold value (either this or "nextremes" must be given but not both) 
}
\item{nextremes}{
the number of upper extremes to be used (either this or "threshold" must be given but not both) 
}
\item{method}{
whether parameters should be estimated by the maximum likelihood method "ml" or 
the probability-weighted moments method "pwm" 
}
\item{information}{
whether standard errors should be calculated with "observed" or "expected" information. 
This only applies to maximum likelihood method; for "pwm" method "expected" information 
is used if possible. 
}
}
\value{
a list containing parameter estimates, standard errors and details of the fit 
}
\details{
see page 278 of QRM; this function uses "nlminb" for ML 
}
\section{References}{
Parameter and quantile estimation for the generalized Pareto distribution, 
JRM Hosking and JR Wallis, Technometrics 29(3), pages 339-349, 1987. 
}
\seealso{
   \code{\link{fit.GPD}},
   \code{\link{fit.GEV}},
   \code{\link{RiskMeasures}}
}
\examples{
data(danish);
losses <- seriesData(danish); 
mod <- fit.GPDb(losses,threshold=10); 
mod;
}
\keyword{methods}

