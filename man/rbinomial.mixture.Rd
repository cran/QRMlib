\name{rbinomial.mixture}
\alias{rbinomial.mixture}
\title{
  Sample Mixed Binomial Distribution 
}
\description{
random generation from mixed binomial distribution 
}
\usage{
rbinomial.mixture(n=1000, m=100, model="probitnorm", ...)
}
\arguments{
\item{n}{
sample size 
}
\item{m}{
vector of numbers of coin flips 
}
\item{model}{
name of mixing distribution: "probitnorm", "logitnorm","beta"," 
}
\item{\dots}{
further parameters of mixing distribution 
}
}
\value{
vector of numbers of successes 
}
\details{
see pages 354-355 and pages 375-377 of QRM 
}
\seealso{
   \code{\link{rbeta}},
   \code{\link{rprobitnorm}},
   \code{\link{rlogitnorm}}
}
\examples{
pi <- 0.04896; #one obligor defaulting pi = .04896 
pi2 <- 0.00321; #two obligors defaulting jointly pi2 = .0031265
beta.pars <- cal.beta(pi,pi2); 
probitnorm.pars <- cal.probitnorm(pi,pi2); 
n <- 1000; 
m <- rep(500,n); 
M.beta <- rbinomial.mixture(n,m,"beta",shape1=beta.pars[1],
         shape2=beta.pars[2]); 
M.probitnorm <- rbinomial.mixture(n,m,"probitnorm",
    mu=probitnorm.pars[1],sigma=probitnorm.pars[2]); 
}
\keyword{distribution}

