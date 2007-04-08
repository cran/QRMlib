\name{rGIG}
\alias{rGIG}
\title{
  Generate Random Vector from Generalized Inverse Gaussian Distribution 
}
\description{
random generation for the generalized inverse Gaussian distribution 
}
\usage{
rGIG(n, lambda, chi, psi, envplot=FALSE, messages=FALSE)
}
\arguments{
\item{n}{
sample size 
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
\item{envplot}{
whether a plot of the rejection envelope should be made; default is FALSE 
}
\item{messages}{
whether a message about rejection rate should be returned; default is FALSE 
}
}
\value{
random sample of size n 
}
\details{
uses a rejection algorithm suggested by Atkinson (1982) 
}
\section{References}{
Atkinson A.C. (1982). The simulation of generalized inverse 
Gaussian and hyperbolic random variables. SIAM Journal on 
Scientific Computing 3(4): 502-515. 
}
\seealso{
   \code{\link{rghyp}},
   \code{\link{rmghyp}}
}
\examples{
\dontrun{
#Create a mean-variance normal mixture of random 
#variables called the generalized hyperbolic
#It is not necessarily elliptical but its univariate 
#version will be.  See p. 78 in QRM.
# This is the GH model.
rghyp <- function(n, lambda, chi, psi, mu=0, gamma=0)
{
  #generate a series of random Generalized Inverse Gaussian 
  #variables: see p. 77 of QRM text
  W <- rGIG(n, lambda, chi, psi);
  # Generate a similar random sequence of standard normals: 
  Z <- rnorm(n);
  #Mix the two distributions using equation 3.25 (p. 77) but 
  #with gamma possibly 0 or a scalar
  sqrt(W) * Z + mu + gamma * W;
}
}
}

\keyword{distribution}

