\name{fit.tcopula.rank}
\alias{fit.tcopula.rank}
\title{
  Fit t Copula Using Rank Correlations 
}
\description{
fits t copula to pseudo-copula data 
}
\usage{
fit.tcopula.rank(Udata, method="Kendall")
}
\arguments{
\item{Udata}{
matrix of pseudo-copula data where rows are vector observations with all values in unit interval 
}
\item{method}{
method to use for calculating rank correlations; default is "Kendall", which is theoretically justified 
}
}
\value{
list containing parameter estimates and details of fit 
}
\details{
see pages 229-231 of QRM 
}
\seealso{
   \code{\link{fit.tcopula}},
   \code{\link{fit.gausscopula}},
   \code{\link{fit.Archcopula2d}}
}

\examples{
\dontrun{
# Multivariate Fitting with Gauss and t: Simulation
# Create an equicorrelation matrix:
P <- equicorr(3,0.6);
set.seed(113);
#Generate a new set of random data from a t-copula (10df) with the same Sigma matrix:
Udatasim2 <- rcopula.t(1000,df=10,Sigma=P);
#Now fit the copula to the simulated data using (Kendall) rank correlations
#and the fit.tcopula.rank() method:
mod.t2 <- fit.tcopula.rank(Udatasim2);
mod.t2;
}
}
\keyword{methods}

