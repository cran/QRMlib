\name{fit.st}
\alias{fit.st}
\title{
  Fit Student t Distribution 
}
\description{
fits univariate Student's t distribution 
}
\usage{
fit.st(data)
}
\arguments{
\item{data}{
vector of data 
}
}
\value{
list containing parameter estimates, standard errors and details of fit 
}
\details{
see page 75 of QRM 
}
\seealso{
   \code{\link{fit.NH}},
   \code{\link{fit.mNH}},
   \code{\link{fit.mst}}
}
\examples{
data(DJ);
Ret.DJ <- mk.returns(DJ); 
window1.start <- timeDate("01/01/1993",format="\%m/\%d/\%Y"); 
window1.end <-  timeDate("12/31/2000",format="\%m/\%d/\%Y"); 
sample1 <- (seriesPositions(Ret.DJ) > window1.start 
     & seriesPositions(Ret.DJ) < window1.end); 
DJ30daily <- Ret.DJ[sample1,]; 
DJ30daily <- 100*seriesData(DJ30daily); 
rseries <- DJ30daily[,"MSFT"]; 
mod.t <- fit.st(rseries); 
\dontrun{
mod.gauss <- fit.norm(rseries);
c(mod.gauss$ll.max, mod.t$ll.max);
mod.NIG <- fit.NH(rseries); 
mod.hyp <- fit.NH(rseries,case="hyp"); 
c(mod.NIG$ll.max, mod.hyp$ll.max); 
}
}
\keyword{methods}

