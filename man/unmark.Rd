\name{unmark}
\alias{unmark}
\title{Unmark Point Process}
\description{
  strips marks from a marked point process   
}
\usage{
unmark(PP)
}
\arguments{
  \item{PP}{a point process object of class PP}
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\seealso{
   \code{\link{fit.sePP}},
   \code{\link{fit.seMPP}},
   \code{\link{extremalPP}}
} 
\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500); 
window <- (seriesPositions(sp500.nreturns) > 
        timeDate("12/31/1995",format = "\%m/\%d/\%Y")); 
sp500.nreturns <- sp500.nreturns[window]; 
tmp <- extremalPP(sp500.nreturns,ne=100); 
tmp2 <- unmark(tmp); 
}
\keyword{methods}

