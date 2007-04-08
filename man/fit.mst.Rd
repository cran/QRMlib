\name{fit.mst}
\alias{fit.mst}
\title{
  Fit Multivariate Student t Distribution 
}
\description{
fits multivariate Student's t distribution using variant of EM algorithm 
}
\usage{
fit.mst(data=data.t.5d, nit=2000, tol=1e-10)
}
\arguments{
\item{data}{
matrix of data where rows are vector observations 
}
\item{nit}{
number of iterations of EM-type algorithm 
}
\item{tol}{
tolerance of improvement for stopping iteration 
}
}
\value{
list containing parameter estimates, standard errors and details of fit 
}
\details{
see page 75 of QRM 
}
\seealso{
   \code{\link{fit.mNH}},
   \code{\link{fit.NH}},
   \code{\link{fit.st}}
}
\examples{
data(DJ);
Ret.DJ <- mk.returns(DJ); 
window1.start <- timeDate("01/01/1993",format="\%m/\%d/\%Y"); 
window1.end <-  timeDate("12/31/2000",format="\%m/\%d/\%Y"); 
sample1 <- (seriesPositions(Ret.DJ) > window1.start 
            & seriesPositions(Ret.DJ) < window1.end); 
selection1 <- c("AXP","EK","BA","C","KO","MSFT",
              "HWP","INTC","JPM","DIS"); 
DJ30daily <- Ret.DJ[sample1,selection1]; 
DJ30weekly <-  aggregateWeeklySeries(DJ30daily, FUNC= colSums); 
mod.t <- fit.mst(DJ30weekly);
\dontrun{
mod.GAUSS <- fit.norm(DJ30weekly);
mod.NIG <- fit.mNH(DJ30weekly,symmetric=FALSE,case="NIG"); 
mod.HYP <- fit.mNH(DJ30weekly,symmetric=FALSE,case="hyp"); 
mod.NIGs <- fit.mNH(DJ30weekly,symmetric=TRUE,case="NIG"); 
mod.HYPs <- fit.mNH(DJ30weekly,symmetric=TRUE,case="hyp"); 
round(c(mod.GAUSS$ll.max,mod.t$ll.max,mod.NIGs$ll.max,
    mod.HYPs$ll.max,mod.NIG$ll.max,mod.HYP$ll.max),1); 
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

