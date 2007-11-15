\name{EMupdate}
\alias{EMupdate}
\title{
  EM Update Step for Generalized Hyperbolic Estimation
}
\description{
updates estimates of location (mu), dispersion (Sigma) and skewness 
(gamma) parameters in EM estimation of multivariate generalized 
hyperbolic distributions 
}
\usage{
EMupdate(data, mix.pars, mu, Sigma, gamma, symmetric, 
        scaling=TRUE, kvalue=1)
}
\arguments{
\item{data}{
data matrix 
}
\item{mix.pars}{
current values of lambda, chi and psi 
}
\item{mu}{
current value of mu 
}
\item{Sigma}{
current value of sigma 
}
\item{gamma}{
current value of gamma 
}
\item{symmetric}{
logical variable for elliptically symmetric case 
}
\item{scaling}{
do we scale determinant of Sigma to be fixed value? 
}
\item{kvalue}{
value of determinant in the case of scaling. If not scaling, 
you do not need to pass this parameter but can let R set its default.
}
}
\value{
a list with updated estimates of mu (location), Sigma (dispersion) and gamma (skewness)
}
\details{
See pp 81-83 of QRM; in that case k is the determinant of the sample covariance matrix.
\sQuote{EM} stands for the \dQuote{Expectation-Maximization} type of algorithm used to
fit proposed multivariate hyperbolic models to actual data.
 
}
\seealso{
   \code{\link{fit.mNH}}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{multivariate}

