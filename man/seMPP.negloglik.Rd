\name{seMPP.negloglik}
\alias{seMPP.negloglik}
\title{Marked Self-Exciting Point Process Log-Likelihood}
\description{
  evaluates negative log-likelihood of a marked self-exciting point process model;
  this will be objective function massed to nlminb() or optim().
}
\usage{
seMPP.negloglik(theta, PP, case, markdens)
}
\arguments{
  \item{theta}{vector of parameters of self-exciting model}
  \item{PP}{point-process data}
  \item{case}{a numerical variable coding whether Hawkes or ETAS forms are used 
             and whether marks may have an influence on future points }
  \item{markdens}{name of density for marks; currently must be "GPD"}
}

\value{value of log-likelihood}

\seealso{
 \code{\link{fit.seMPP}},
 \code{\link{fit.sePP}}
}
\examples{
\dontrun{
#Example of using seMPP.negloglik as objective function passed 
#to optimizer function
fit.seMPP <- function(PP,markdens="GPD",model="Hawkes",
       mark.influence=TRUE,predictable=FALSE,std.errs=FALSE)
{
  if (class(PP) != "MPP") stop("Not marked point process data");
  marks <- PP$marks;
  groundmod <- fit.sePP(PP,model,mark.influence=TRUE,std.errs=FALSE);
  #lines removed here...
  if(predicatable)
  {
    theta <- c(par.ests,0); 
    fit <- nlminb(start=theta, objective=seMPP.negloglik, 
            PP=PP,case=case, markdens=markdens);
    #Lines removed here ...
  }
}
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{optimize}

