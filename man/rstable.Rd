\name{rstable}
\alias{rstable}
\title{
  Stable Distribution 
}
\description{
random sample from stable distribution 
}
\usage{
rstable(n, alpha, beta=1)
}
\arguments{
\item{n}{
sample size 
}
\item{alpha}{
scalar parameter strictly larger than 0 and smaller 
than 2 (but avoid alpha=1) 
}
\item{beta}{
scalar parameter between -1 and 1 
}
}
\value{
sample of size n 
}
\details{
see pages 224 and 498 of QRM; default value beta=1 combined 
with an alpha value less than 1 gives positive stable distribution 
which we require for Gumbel copula generation; the case alpha=1 
has not been implemented 
}
\section{References}{
Forthcoming John Nolan Book; see Definition 1.8 and Theorem 1.19 
}
\seealso{
   \code{\link{rcopula.gumbel}}
}
\examples{
\dontrun{
#Use rstable() method in copula simulation function
rcopula.Gumbel2Gp <- function(n = 1000, 
       gpsizes =c(2,2), theta =c(2,3,5))
{
  Y <- rstable(n,1/theta[1])*(cos(pi/(2*theta[1])))^theta[1];
  innerU1 <- rcopula.gumbel(n,theta[2]/theta[1],gpsizes[1]);
  innerU2 <- rcopula.gumbel(n,theta[3]/theta[1],gpsizes[2]);
  U <- cbind(innerU1,innerU2);
  Y <- matrix(Y, nrow = n, ncol = sum(gpsizes));                               
  out <- exp( - ( - log(U)/Y)^(1/theta[1]));
  return(out);
}
}

}
\keyword{distribution}

