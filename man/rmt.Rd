\name{rmt}
\alias{rmt}
\title{
  Multivariate t 
}
\description{
generates random sample from multivariate t 
}
\usage{
rmt(n, df=4, Sigma=equicorr(d, rho), mu=rep(0, d), d=2, rho=0.7)
}
\arguments{
\item{n}{
number of realizations 
}
\item{df}{
degrees of freedom 
}
\item{Sigma}{
a dispersion matrix 
}
\item{mu}{
a mean vector 
}
\item{d}{
dimension of distribution 
}
\item{rho}{
correlation value to build equicorrelation matrix 
}
}
\value{
an n by d matrix 
}
\details{
function is set up to quickly simulate equicorrelation 
structures by specifying d and rho 
}
\seealso{
   \code{\link{rmnorm}},
   \code{\link{equicorr}},
   \code{\link{rmghyp}}
}
\examples{
\dontrun{
ndata <- rmnorm(1000,df=4,rho=0.7,d=3); 
Sigma <- diag(c(3,4,5)) \%*\% equicorr(3,0.6) 
         \%*\% diag(c(3,4,5)); 
mu <- c(1,2,3); 
tdata <- rmt(1000,4,Sigma,mu); 
mod1 <- fit.mst(tdata);
}
}
\keyword{distribution}

