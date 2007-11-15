\name{plotFittedGPDvsEmpiricalExcesses}
\alias{plotFittedGPDvsEmpiricalExcesses}
\title{Graphically Compare Empirical Distribution of Excesses and GPD Fit}
\description{
  Build a graph which plots the GPD fit of excesses over a threshold u
and the corresponding empirical distribution function for observed excesses.
}
\usage{
plotFittedGPDvsEmpiricalExcesses(data, threshold = NA, nextremes = NA)
}
\arguments{
  \item{data}{data vector or times series; to be safe, pass data slot of timeSeries}
  \item{threshold}{a threshold value (either this or "nextremes" must be given but not both)}
  \item{nextremes}{the number of upper extremes to be used (either this or "threshold" must be given but not both)}
}
\details{
  See graphs 7.4(c) and 7.5(c) in QRM, pp. 281-2. 
}
\value{
  a plot showing empirical cdf of excesses vs points fitted to the estimated GPD for excesses
}

\seealso{
   \code{\link{fit.GPD}},
   \code{\link{plotTail}},
   \code{\link{MEplot}},
   \code{\link{xiplot}}
}

\examples{
data(danish);
plotFittedGPDvsEmpiricalExcesses(danish@Data, nextremes=109);
plotFittedGPDvsEmpiricalExcesses(danish@Data,threshold=10);
}
\keyword{hplot}

