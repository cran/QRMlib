\name{dmnorm}
\alias{dmnorm}
\title{
  Multivariate Normal Density 
}
\description{
evaluates multivariate normal density 
}
\usage{
dmnorm(x, mu, Sigma, logvalue=FALSE)
}
\arguments{
\item{x}{
matrix with n rows and d columns; density is evaluated at each vector of row values 
}
\item{mu}{
mean vector 
}
\item{Sigma}{
covariance matrix 
}
\item{logvalue}{
should log density be returned; default is FALSE 
}
}
\value{
vector of length n containing values of density or log-density 
}
\seealso{
   \code{\link{dmt}},
   \code{\link{dmghyp}}
}
\examples{
### Normal distribution: visualization, simulation, estimation
BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7));
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{distribution}

