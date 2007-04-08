\name{volfunction}
\alias{volfunction}
\title{Self-Excitement Function}
\description{
  calculates a self-excitement function for use in the negloglik() methods
 used in fit.sePP() and fit.seMPP()
}
\usage{
volfunction(anytimes, times, marks, theta, model)
}
\arguments{
  \item{anytimes}{vector of times at which to calculate self-excitement function}
  \item{times}{times of point events}
  \item{marks}{marks associated with point events}
  \item{theta}{parameters of self-excitement function}
  \item{model}{model number}
}
\details{
  see page 306 of QRM  
}
\value{
  a vector of same length as "anytimes"
}

\seealso{
   \code{\link{fit.sePP}},
   \code{\link{fit.seMPP}}
}
\examples{
\dontrun{
seMPP.negloglik <- function(theta, PP, case, markdens) 
 {
   theta <- abs(theta);
   times <- PP$times;
   marks <- PP$marks;
   endtime <- PP$endtime;
   starttime <- PP$starttime;
   mu <- theta[1];
   phi <- theta[2];
   voltheta <- theta[-c(1,2,(length(theta)-2),
        (length(theta)-1),length(theta))];
   evol <- volfunction(times,times,marks,voltheta,case);
   lambda.contrib <- mu+phi*evol
   # ... remaining lines ommitted here
 }
}
}
\keyword{methods}

