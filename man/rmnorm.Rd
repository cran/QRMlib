\name{rmnorm}
\alias{rmnorm}
\title{
  Multivariate Normal Random Sample
}
\description{
generates random sample from multivariate normal 
}
\usage{
rmnorm(n, Sigma=equicorr(d, rho), mu=rep(0, d), d=2, rho=0.7)
}
\arguments{
\item{n}{
number of realizations 
}
\item{Sigma}{
a covariance matrix 
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
function is set up to quickly simulate equicorrelation structures 
by specifying d and rho 
}
\seealso{
   \code{\link{rmt}},
   \code{\link{equicorr}}
}
\examples{
ndata <- rmnorm(1000,rho=0.7,d=3); 
\dontrun{
Sigma <- diag(c(3,4,5)) \%*\% equicorr(3,0.6) \%*\% diag(c(3,4,5)); 
mu <- c(1,2,3); 
ndata <- rmnorm(1000,Sigma,mu);
fit.norm(ndata);
} 
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

