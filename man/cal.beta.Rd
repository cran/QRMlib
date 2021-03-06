\name{cal.beta}
\alias{cal.beta}
\title{
  Calibrate Beta Mixture of Bernoullis 
}
\description{
calibrates a beta mixture distribution on unit interval to give an 
exchangeable Bernoulli mixture model with prescribed default and 
joint default probabilities 
}
\usage{
cal.beta(pi1=0.1837, pi2=0.0413)
}
\arguments{
\item{pi1}{
default probability 
}
\item{pi2}{
joint default probability 
}
}
\value{
parameters a and b of beta mixing distribution 
}
\details{
see pages 354-355 in QRM 
}
\seealso{
   \code{\link{cal.claytonmix}},
   \code{\link{cal.probitnorm}},
   \code{\link{rbinomial.mixture}}
}

\examples{
pi.B <- 0.2; pi2.B <- 0.05 
probitnorm.pars <- cal.probitnorm(pi.B,pi2.B) 
probitnorm.pars 
beta.pars <- cal.beta(pi.B,pi2.B) 
beta.pars 
claytonmix.pars <- cal.claytonmix(pi.B,pi2.B) 
claytonmix.pars 
q <- (1:1000)/1001; 
q <- q[q<0.25]; 
p.probitnorm <- pprobitnorm(q,probitnorm.pars[1],
               probitnorm.pars[2]); 
p.beta <- pbeta(q, beta.pars[1], beta.pars[2]); 
p.claytonmix <- pclaytonmix(q,claytonmix.pars[1],
                  claytonmix.pars[2]); 
scale <- range((1-p.probitnorm),(1-p.beta),(1-p.claytonmix)); 
plot(q, (1 - p.probitnorm), type = "l", log = "y", xlab = "q", 
           ylab = "P(Q>q)",ylim=scale); 
lines(q, (1 - p.beta), col = 2); 
lines(q, (1 - p.claytonmix), col = 3); 
legend("topright", c("Probit-normal", "Beta", "Clayton-Mixture"), 
          lty=rep(1,3),col = (1:3)) 
}
\keyword{models}

