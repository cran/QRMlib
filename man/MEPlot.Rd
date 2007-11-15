\name{MEplot}
\alias{MEplot}
\title{
  Sample Mean Excess Plot 
}
\description{
Plots sample mean excesses over increasing thresholds.  
}
\usage{
MEplot(data, omit = 3, labels=TRUE, ...)
}
\arguments{
\item{data}{
data vector or Data slot from a time series (e.g tS@Data); do not pass entire timeSeries
}
\item{omit}{
number of upper plotting points to be omitted  
}
\item{labels}{
whether or not axes are to be labelled 
}
\item{\dots}{
further parameters of MEplot function 
}
}
\details{
An upward trend in plot shows heavy-tailed behaviour. In particular, a 
straight line with positive gradient above some threshold is a sign 
of Pareto behaviour in tail. A downward trend shows thin-tailed 
behaviour whereas a line with zero gradient shows an exponential 
tail. 
Because upper plotting points are the average of a handful of extreme 
excesses, these may be omitted for a prettier plot. 
}
\seealso{
   \code{\link{fit.GPD}}
}
\examples{
# Sample mean excess plot of heavy-tailed Danish fire insurance data
data(danish);
MEplot(data=danish@Data); 
}
\keyword{hplot}
