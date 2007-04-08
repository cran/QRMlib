\name{momest}
\alias{momest}
\title{
  Moment Estimator of Default Probabilities 
}
\description{
calculates moment estimator of default probabilities and joint default probabilities 
for a homogeneous group  
}
\usage{
momest(data, trials, limit=10.)
}
\arguments{
\item{data}{
vector of numbers of defaults in each time period 
}
\item{trials}{
vector of group sizes in each time period 
}
\item{limit}{
maximum order of joint default probability to estimate 
}
}
\value{
vector of default probability and joint default probability estimates 
}
\details{
first returned value is default probability estimate; second value is estimate of joint default probability for two firms; and so on. See pages 375-376 in QRM 
}
\seealso{
   \code{\link{fit.binomialBeta}},
   \code{\link{fit.binomialLogitnorm}},
   \code{\link{fit.binomialProbitnorm}}
}
\examples{
\dontrun{
#MODEL RISK See especially Section 8.4.6 on p. 364 of QRM book
spdata.raw;
attach(spdata.raw);
#momest() is an internal function in functionsCredit.R to 
#calculate moment estimators for default probabilities. The first 
#parameter input is a vector containing the number of defaults in 
#each time period; the 2nd parameter input is a vector containing the
#number of credits in the group during the time period.
momest(Bdefaults,Bobligors);
#The values calculated from momest(Bdefaults,Bobligors) are the 
#parameter estimates shown in Table 8.6, p.365 of QRM book under the 
#model column labeled 'B'
#The first value returned is the probability of a single default.
pi.B <- momest(Bdefaults, Bobligors)[1];  #one obligor defaulting pi = .04896
#second value returned is probability of joint default probability for two firms.
pi2.B <- momest(Bdefaults, Bobligors)[2]; #two obligors defaulting jointly pi2 = .0031265
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{htest}

