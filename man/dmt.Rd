\name{dmt}
\alias{dmt}
\title{
  Multivariate Student t Density 
}
\description{
evaluates multivariate Student t density 
}
\usage{
dmt(x, nu, mu, Sigma, logvalue=FALSE)
}
\arguments{
\item{x}{
matrix with n rows and d columns; density is evaluated at each vector of row values 
}
\item{nu}{
degree of freedom parameter 
}
\item{mu}{
location vector 
}
\item{Sigma}{
dispersion matrix 
}
\item{logvalue}{
should log density be returned; default is FALSE 
}
}
\value{
vector of length n containing values of density or log-density 
}
\seealso{
   \code{\link{dmnorm}},
   \code{\link{dmghyp}}
}
\examples{
### t distribution: visualization, simulation, estimation
BiDensPlot(func=dmt,xpts=c(-4,4),ypts=c(-4,4),mu=c(0,0),
            Sigma=equicorr(2,-0.7),nu=4);
}
\keyword{distribution}

