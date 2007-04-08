\name{qst}
\alias{qst}
\title{
  Student's t Distribution (3 parameter) 
}
\description{
Quantiles for 3-parameter version of Student's t distribution 
}
\usage{
qst(p, df, mu=0, sigma=1, scale=FALSE)
}
\arguments{
\item{p}{
vector of probabilities 
}
\item{df}{
vector of degrees of freedom 
}
\item{mu}{
vector of location parameters 
}
\item{sigma}{
vector of scale parameters 
}
\item{scale}{
whether distribution should be scaled so that mu and sigma 
are mean and standard deviation; default is FALSE 
}
}
\value{
quantiles for Student's t distribution 
}
\seealso{
   \code{\link{qt}},
   \code{\link{ESst}}
}
\examples{
#Set up the quantile probabilities
p <- c(0.90,0.95,0.975,0.99,0.995,0.999,0.9999,0.99999,0.999999);
sigma <- 0.2*10000/sqrt(250);
#Now look at VaR for student t with 4 degrees of freedom:
VaR.t4 <- qst(p,4,sigma=sigma,scale=TRUE);
}
\keyword{distribution}

