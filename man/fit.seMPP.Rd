\name{fit.seMPP}
\alias{fit.seMPP}
\title{Fit Marked Self-Exciting Point Process}
\description{
  fits marked self-exciting process to a point process object of class MPP  
}
\usage{
fit.seMPP(PP, markdens = "GPD", model = "Hawkes", mark.influence = TRUE, 
predictable = FALSE, std.errs = FALSE)
}
\arguments{
  \item{PP}{a point process object of class MPP  }
  \item{markdens}{ name of density of mark distribution; currently must be "GPD" }
  \item{model}{ name of self-exciting model: Hawkes or ETAS }
  \item{mark.influence}{whether marks of marked point process may influence the self-excitement}
  \item{predictable}{whether previous events may influence the scaling of mark distribution}
  \item{std.errs}{whether standard errors should be computed VALUE}
}
\details{
  see pages 307-309 of QRM 
}
\value{
  a fitted self-exciting process object of class sePP 
}

\seealso{
   \code{\link{fit.sePP}},
   \code{\link{plot.sePP}},
   \code{\link{stationary.sePP}}
}
\examples{
data(sp500);
sp500.nreturns <- -mk.returns(sp500); 
window <- (seriesPositions(sp500.nreturns) > 
            timeDate("12/31/1995",format = "\%m/\%d/\%Y")); 
sp500.nreturns <- sp500.nreturns[window]; 
tmp <- extremalPP(sp500.nreturns,ne=100); 
mod3a <- fit.seMPP(tmp,mark.influence=FALSE,std.errs=TRUE);
\dontrun{ 
mod3b <- fit.seMPP(tmp,mark.influence=TRUE,std.errs=TRUE); 
mod3c <- fit.seMPP(tmp,model="ETAS",mark.influence=FALSE,std.errs=TRUE); 
mod3d <- fit.seMPP(tmp,model="ETAS",mark.influence=TRUE,std.errs=TRUE); 
mod4a <- fit.seMPP(tmp,mark.influence=FALSE,predictable=TRUE,
           std.errs=TRUE); 
mod4b <- fit.seMPP(tmp,mark.influence=TRUE,predictable=TRUE,
         std.errs=TRUE); 
mod4c <- fit.seMPP(tmp,model="ETAS",mark.influence=FALSE,
                   predictable=TRUE,std.errs=TRUE); 
mod4d <- fit.seMPP(tmp,model="ETAS",mark.influence=TRUE,
                   predictable=TRUE,std.errs=TRUE);
}
}
\keyword{methods }
