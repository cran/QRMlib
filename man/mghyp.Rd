\name{mghyp}
\alias{mghyp}
\alias{dmghyp}
\alias{rmghyp}
\title{
  Multivariate Generalized Hyperbolic Distribution 
}
\description{
Density and random number generation for density of multivariate 
generalized hyperbolic distribution in new QRM (Chi-Psi-Sigma- Gamma) 
parameterization. Note Sigma is the dispersion matrix. See pp. 77-81.
}
\details{
See page 78 in QRM for joint density formula (3.30) with Sigma a d-dimensional
dispersion matrix (d > 1) consistent with a multivariate distribution). This is 
a more intuitive parameterization of the alpha-beta-delta model used by Blaesild 
(1981) in earlier literature since it associates all parameters with mixtures of
both mean and variance.  Here gamma is assumed equal to 0 so we have a normal 
variance mixture where the mixing variable W has a GIG generalized inverse gaussian) 
distribution with parameters lambda, chi, psi. This thickens the tail.

If gamma exceeds zero, we have a normal mean-variance mixture where the mean 
is also perturbed to equal mu + (W * gamma) which introduces ASYMMETRY as well. 

The default d=2 for the random generator gives a two-dimensional matrix of n values.

See pp. 77-81 of QRM and appendix A.2.5 for details.
}
\usage{
dmghyp(x, lambda, chi, psi, mu, Sigma, gamma, logvalue=FALSE)
rmghyp(n, lambda, chi, psi, Sigma=equicorr(d, rho), mu=rep(0, d), 
     gamma=rep(0, d), d=2, rho=0.7)
}
\arguments{
\item{x}{
matrix with n rows and d columns; density is evaluated at each vector of row values 
}
\item{lambda}{
scalar parameter 
}
\item{chi}{
scalar parameter 
}
\item{psi}{
scalar parameter 
}
\item{mu}{
location vector 
}
\item{Sigma}{
dispersion matrix 
}
\item{d}{
dimension of distribution 
}
\item{rho}{
correlation value to build equicorrelation matrix 
}
\item{gamma}{
vector of skew parameters 
}
\item{logvalue}{
should log density be returned; default is FALSE 
}
\item{n}{length of vector}
}
\value{values of density or log-density or randomly generated values}
\note{
See page 78 in QRM; if gamma is a zero vector distribution is elliptical 
and dsmghyp is called.
If lambda = (d+1)/2, we drop generalized and call the density a d-dimensional
hyperbolic density.
If lambda = 1, the univariate marginals are one-dimensional hyperbolics.
If lambda = -1/2, distribution is NIG (normal inverse gaussian).
If lambda greater than 0 and chi = 0, we get the VG (variance gamma)
If we can define a constant nu such that
  lambda = (-1/2)*nu AND chi = nu
then we have a multivariate skewed-t distribution.
See p. 80 of QRM for details. 
}
\seealso{
   \code{\link{dsmghyp}},
   \code{\link{dmt}},
   \code{\link{dmnorm}}
}
\examples{
Sigma <- diag(c(3,4,5)) \%*\% equicorr(3,0.6) \%*\% diag(c(3,4,5)); 
mu <- c(1,2,3); 
ghdata <- rmghyp(n=1000,lambda=0.5,chi=1,psi=1,Sigma,mu); 
### (Multivariate generalized) Hyperbolic distribution: visualization with 
# PERSPECTIVE or CONTOUR plots
par(mfrow=c(2,2));
ll <- c(-4,4);
#pass the multivariate generalized hyperbolic density to be plotted:
BiDensPlot(func=dmghyp,xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),
             lambda=1,chi=1,psi=1,gamma=c(0,0));
BiDensPlot(func=dmghyp,type="contour",xpts=ll,ypts=ll,mu=c(0,0),
           Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0,0));
BiDensPlot(func=dmghyp,xpts=ll,ypts=ll,mu=c(0,0),
           Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0.5,-0.5));
BiDensPlot(func=dmghyp,type="contour",xpts=ll,ypts=ll,mu=c(0,0),
           Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0.5,-0.5));
par(mfrow=c(1,1));
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

