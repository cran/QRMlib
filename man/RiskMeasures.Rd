\name{RiskMeasures}
\alias{RiskMeasures}
\title{
  Calculate Risk Measures from GPD Fit 
}
\description{
calculates risk measures like VaR and expected shortfall 
based on a generalized Pareto model fitted to losses over 
a high threshold 
}
\usage{
RiskMeasures(out, p)
}
\arguments{
\item{out}{
results of a GPD fit to excesses over high thresholds 
}
\item{p}{
vector of probability levels for risk measures 
}
}
\value{
matrix with quantile and shortfall estimates for each 
probability level 
}
\details{
see pages 282-284 of QRM 
}
\seealso{
   \code{\link{fit.GPD}},
   \code{\link{showRM}}
}
\examples{
data(danish);
out <- fit.GPD(danish,10); 
RiskMeasures(out,c(0.99,0.999)); 
}
\keyword{misc}

