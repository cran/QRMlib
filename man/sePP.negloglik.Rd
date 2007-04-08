\name{sePP.negloglik}
\alias{sePP.negloglik}
\title{ Self-Exciting Point Process Log-Likelihood}
\description{
  evaluates negative log-likelihood of a self-exciting point 
  process model (unmarked) 
}
\usage{
sePP.negloglik(theta, PP, case)
}
\arguments{
  \item{theta}{parameters of self-exciting model}
  \item{PP}{point-process data}
  \item{case}{a numerical variable coding whether Hawkes 
   or ETAS forms are used and whether marks may have an 
   influence on future points}
}
\value{
  value of log-likelihood
}

\seealso{
   \code{\link{fit.sePP}},
   \code{\link{fit.seMPP}}
}
\examples{
\dontrun{
#Example of using sePP.negloglik as objective function passed 
#to optimizer function
fit.sePP <- function(PP,markdens="GPD",model="Hawkes",
      mark.influence=T,predictable=F,std.errs=F)
{
  #lines removed here...
  fit <- nlminb(start=theta, objective=sePP.negloglik, PP=PP,case=case);
  par.ests <- fit$par;
  par.ests <- abs(par.ests)
  ll.max <- -sePP.negloglik(par.ests,PP,case)
  #Lines removed here ...
}
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{optimize}

