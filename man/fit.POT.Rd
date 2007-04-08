\name{fit.POT}
\alias{fit.POT}
\title{ Peaks-over-Threshold Model}
\description{
fits the POT (peaks-over-threshold) model to a point process 
of class PP or MPP 
}
\usage{
fit.POT(PP, markdens = "GPD")
}
\arguments{
  \item{PP}{a point process object of class PP or MPP}
  \item{markdens}{(optional) 
   name of density of mark distribution, currently must be "GPD"}
}
\details{
  see pages 301-305 of QRM.  Note that if point process is of 
  class PP then function simply esitmates the rate of a homogeneous 
  Poisson process.  
}
\value{
  a list containing parameters of fitted POT model  
  \item{par.ests}{vector of parameter estimates}
  \item{par.ses}{vector of parameter std deviations}
  \item{ll.max}{logliklihood maximum}
}
\references{see pages 301-305 of QRM }

}
\seealso{
   \code{\link{fit.GPD}},
   \code{\link{extremalPP}}
}
\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500); 
window <- (seriesPositions(sp500.nreturns) > 
         timeDate("12/31/1995",format="\%m/\%d/\%Y")); 
sp500.nreturns <- sp500.nreturns[window];
tmp <- extremalPP(sp500.nreturns,ne=100); 
mod1 <- fit.POT(tmp); 
\dontrun{
mod1b <- fit.POT(unmark(tmp)); 
}
}
\keyword{methods}

