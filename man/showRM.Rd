\name{showRM}
\alias{showRM}
\title{
  Show Risk Measure Estimates on Tailplot 
}
\description{
shows estimates of risk measures (like VaR and ES) on a tailplot 
}
\usage{
showRM(object, alpha, RM="VaR", extend=2, ci.p=0.95, like.num=50.)
}
\arguments{
\item{object}{
results of fit.GPD 
}
\item{alpha}{
probability level 
}
\item{RM}{
risk measure, VaR or ES 
}
\item{extend}{
how far to extend picture; x-axis extends to this value times the largest observation 
}
\item{ci.p}{
confidence level for confidence interval 
}
\item{like.num}{
number of evaluations of profile likelihood 
}
}
\value{
point estimate and confidence interval for risk measure 
}
\section{Side Effects}{
plotTail is called 
}
\details{
see pages 282-284 in QRM 
}
\seealso{
   \code{\link{plotTail}},
   \code{\link{fit.GPD}},
   \code{\link{RiskMeasures}}
}
\examples{
\dontrun{
data(danish);
#Fit the GPD using MLE a
mod <- fit.GPD(danish,10); 
showRM(mod,0.99,RM="VaR");
showRM(mod,0.99,RM="ES");
showRM(mod,0.995,RM="VaR");
showRM(mod,0.995,RM="ES");
}
}
\keyword{hplot}

