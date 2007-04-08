\name{fit.NH}
\alias{fit.NH}
\title{
  Fit NIG or Hyperbolic Distribution 
}
\description{
fits univariate NIG or hyperbolic distribution 
}
\usage{
fit.NH(data, case="NIG", symmetric=FALSE, se=FALSE)
}
\arguments{
\item{data}{
vector of data 
}
\item{case}{
whether NIG  ("NIG") or hyperbolic ("hyp"); default is NIG 
}
\item{symmetric}{
whether symmetric or asymmetric; default is FALSE 
}
\item{se}{
whether standard errors should be calculated 
}
}
\value{
list containing parameter estimates, standard errors and details of fit 
}
\details{
See pages 78-80 of QRM. Case \sQuote{NIG} sets lambda to -1/2; case \sQuote{hyp}
sets lambda to 1; no other cases are allowed.
}
\seealso{
   \code{\link{fit.st}},
   \code{\link{fit.mNH}},
   \code{\link{fit.mst}}
}
\examples{
data(DJ);
Ret.DJ <- mk.returns(DJ); 
window1.start <- timeDate("01/01/1993", format="\%m/\%d/\%Y"); 
window1.end <-  timeDate("12/31/2000",format="\%m/\%d/\%Y"); 
sample1 <- (seriesPositions(Ret.DJ) > window1.start 
      & seriesPositions(Ret.DJ) < window1.end); 
DJ30daily <- Ret.DJ[sample1,]; 
DJ30daily <- 100*seriesData(DJ30daily); 
rseries <- DJ30daily[,"MSFT"]; 
mod.NIG <- fit.NH(rseries);
\dontrun{
mod.gauss <- fit.norm(rseries);
mod.t <- fit.st(rseries); 
mod.hyp <- fit.NH(rseries,case="hyp"); 
c(mod.gauss$ll.max,mod.t$ll.max,mod.NIG$ll.max,mod.hyp$ll.max); 
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

