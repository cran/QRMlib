\name{ESst}
\alias{ESst}
\title{
  Expected Shortfall for Student t Distribution 
}
\description{
calculates expected shortfall for Student t distribution 
}
\usage{
ESst(p, df, mu=0, sigma=1, scale=FALSE)
}
\arguments{
\item{p}{
probability level 
}
\item{df}{
degrees of freedom 
}
\item{mu}{
mean 
}
\item{sigma}{
standard deviation 
}
\item{scale}{
should t distribution be scaled to have variance one? 
}
}
\value{
expected shortfall 
}
\details{
see page 45 of QRM 
}
\seealso{
   \code{\link{ESnorm}}
}
\examples{
#Set up the quantile probabilities
p <- c(0.90,0.95,0.975,0.99,0.995,0.999,0.9999,0.99999,0.999999);
sigma <- 0.2*10000/sqrt(250);
#Now look at Expected Shortfall for student t with 4 degrees of freedom:
ES.t4 <- ESst(p,4,sigma=sigma,scale=TRUE);
ESst(c(0.95,0.99),4);
}
\keyword{distribution}

