\name{fit.sePP}
\alias{fit.sePP}
\title{Fit Self-Exciting Process}
\description{
  fits fits self-exciting process to a point process object of class 
  PP (unmarked) or MPP (marked) 
}
\usage{
fit.sePP(PP, model = "Hawkes", mark.influence = TRUE, std.errs = FALSE)
}
\arguments{
  \item{PP}{a point process object of class PP (unmarked) or MPP (marked)}
  \item{model}{(optional)name of self-exciting model: Hawkes or ETAS}
  \item{mark.influence}{(optional)whether marks of marked point process 
    may influence the self-excitement }
  \item{std.errs}{(optional) whether standard errors should be computed VALUE}
}
\details{
  see pages 306-307 of QRM 
}
\value{
  a fitted self-exciting process object of class sePP 
}

\seealso{
   \code{\link{fit.seMPP}},
   \code{\link{plot.sePP}},
   \code{\link{stationary.sePP}}
}
\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500); 
window <- (seriesPositions(sp500.nreturns) > 
    timeDate("12/31/1995",format="\%m/\%d/\%Y")); 
sp500.nreturns <- sp500.nreturns[window]; 
tmp <- extremalPP(sp500.nreturns,ne=100); 
mod2a <- fit.sePP(tmp,mark.influence=FALSE,std.errs=TRUE); 
}
\keyword{methods}
