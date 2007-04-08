\name{ghyp}
\alias{dghyp}
\alias{rghyp}
\title{
  Univariate Generalized Hyperbolic Distribution 
}
\description{
Density and random number generation for univariate generalized hyperbolic 
distribution in new QRM (Chi-Psi-Gamma) parameterization. (The dispersion matrix
Sigma is identically 1, i.e. a scalar 1.) See pp. 77-81 in QRM.
}
\details{
See page 78 in QRM for joint density formula (3.30) with Sigma (dispersion matrix)
the identity and d=1 (meaning a univariate distribution) applies. 

The univariate QRM parameterization is defined in terms of parameters chi-psi-gamma
instead of the alpha-beta-delta model used by Blaesild (1981) in earlier
literature.  If gamma is 0, we have a normal variance mixture where the mixing 
variable W has a GIG generalized inverse gaussian) distribution with parameters lambda, 
chi, psi. This thickens the tail.

If gamma exceeds zero, we have a normal mean-variance mixture where the mean 
is also perturbed to equal mu + (W * gamma) which introduces ASYMMETRY as well. 

Values for lambda and mu are identical in both QRM and B parameterizations.
Sigma does not appear in the parameter list since in the univariate case its value is 
identically 1.
}
\note{
Density values from dgyhp() should be identical to those from dghypB() 
if the alpha-beta-delta parameters of the B type are translated to the
corresponding gamma-chi-psi parameters of the QRM type by formulas on pp 79-80.
}
\usage{
dghyp(x, lambda, chi, psi, mu=0, gamma=0, logvalue=FALSE) 
rghyp(n, lambda, chi, psi, mu=0, gamma=0)
}
\arguments{
\item{x}{
vector of values at which to evaluate density 
}
\item{n}{
sample size 
}
\item{lambda}{
scalar mixing parameter 
}
\item{chi}{
scalar mixing parameter 
}
\item{psi}{
scalar mixing parameter 
}
\item{mu}{
location parameter 
}
\item{gamma}{
skewness parameter 
}
\item{logvalue}{
should log density be returned; default is FALSE 
}
}
\value{
values of density or log-density (dghyp) or random sample (rghyp) 
}

\seealso{
   \code{\link{dghypB}},
   \code{\link{besselM3}},
   \code{\link{dmghyp}}
}
\examples{
data(DJ);
#Make returns from timeSeries (the default is log-returns). 
#Ret.DJ is a timeSeries class.
Ret.DJ <- mk.returns(DJ);
#In R, cut() method works by selecting data only between the 
#'to' and 'from' dates. Hence  we will use the remaining (cut) data from 
#1993-01-01 to 2000-12-31. MUST use PRIOR DAY on 'from' 
DJ30dailyTS <- cut(Ret.DJ, from="1992-12-31", to="2000-12-31");
DJ30daily <- 100*seriesData(DJ30dailyTS);
#Extract only the Microsoft returns as 'rseries'; remember this is a vector--not a timeSeries
rseries <- DJ30daily[,"MSFT"];
#The default case for fit.NH(() is NIG requiring lambda = -1/2.
mod.NIG <- fit.NH(rseries);
xvals <- seq(from=min(rseries),to=max(rseries),length=100);
yvals.NIG <- dghyp(xvals,lambda=-1/2,chi=mod.NIG$par.ests[1],
  psi=mod.NIG$par.ests[2],mu=mod.NIG$par.ests[3],gamma=mod.NIG$par.ests[4]);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

