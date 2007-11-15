\name{fit.mst}
\alias{fit.mst}
\title{
  Fit Multivariate Student t Distribution 
}
\description{
fits multivariate Student's t distribution using variant of EM algorithm 
}
\usage{
fit.mst(data, nit=2000, tol=1e-10)
}
\arguments{
\item{data}{
matrix of data where rows are vector observations.  A good choice
might be data.t.5d =
rmghyp(n,lambda=(-nu/2),chi=nu,psi=0,Sigma=P,mu=mu,gamma=gamma)
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
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

