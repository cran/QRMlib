\name{fit.mNH}
\alias{fit.mNH}
\title{
  Fit Multivariate NIG or Hyperbolic Distribution 
}
\description{
fits multivariate NIG or hyperbolic distribution using variant of EM algorithm 
}
\usage{
fit.mNH(data, symmetric=FALSE, case="NIG", 
     kvalue=NA, nit=2000, tol=1e-10)
}
\arguments{
\item{data}{
matrix of data where rows are vector observations; common example is data.hyp.5d 
}
\item{symmetric}{
whether symmetric case should be fitted; default is FALSE 
}
\item{case}{
whether NIG ("NIG") or hyperbolic ("hyp") should be fitted 
}
\item{kvalue}{
value to which to constrain determinant of dispersion matrix 
}
\item{nit}{
maximum number of iterations 
}
\item{tol}{
tolerance for convergence 
}
}
\value{
list containing parameter estimates, standard errors and details of fit 
}
\details{
see pages 81-83 in QRM 
}
\seealso{
   \code{\link{fit.mst}},
   \code{\link{fit.NH}},
   \code{\link{EMupdate}},
   \code{\link{MCECMupdate}},
   \code{\link{MCECM.Qfunc}}
}
\examples{
data(DJ);
Ret.DJ <- mk.returns(DJ); 
window1.start <- timeDate("01/01/1993",format="\%m/\%d/\%Y"); 
window1.end <-  timeDate("12/31/2000", format="\%m/\%d/\%Y"); 
sample1 <- (seriesPositions(Ret.DJ) > window1.start & 
            seriesPositions(Ret.DJ) < window1.end); 
selection1 <- c("AXP","EK","BA","C","KO","MSFT",
                "HWP","INTC","JPM","DIS"); 
DJ30daily <- Ret.DJ[sample1,selection1]; 
DJ30weekly <- aggregateWeeklySeries(DJ30daily, FUNC= colSums);
mod.NIG <- fit.mNH(DJ30weekly,symmetric=FALSE,case="NIG"); 
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

