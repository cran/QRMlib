\name{QQplot}
\alias{QQplot}
\title{
  Generic Quantile-Quantile Plot 
}
\description{
constructs a quantile-quantile plot against a given reference distribution
(only normal, exponential, and student-t are currently supported) 
}
\usage{
QQplot(data, position=0.5, reference="normal", ...)
}
\arguments{
\item{data}{
vector of data 
}
\item{position}{
determines the plotting positions (see ppoints in R-help) 
}
\item{reference}{
name of reference distribution (only normal, exp, student-t currently allowed) 
}
\item{\dots}{
Any further parameters required by quantile function of reference distribution. 
For example, if reference="exp", you may want to pass 'rate'.  If using the normal
distribution, you may want to pass mu and sigma. If using the student-t, you may
wish to pass the df parameter.
}
}
\details{
  Special forms like ParetoQQ plots can also be created via this function. E.g.,
  to create a ParetoQQ plot, merely pass log(data) in place of data as the first
  parameter and use reference="exp" as the reference distribution.  The ParetoQQ plot
  should provide a linear graph when a log transform of the data is plotted against
  the exponential distribution.  See Beirlant et al, "Statistics of Extremes", Chapter 1.2.1
  for descriptions of various QQ plots.
}
\value{
NULL returned 
}
\section{Side Effects}{
QQ-plot is created on graphical device 
}
\seealso{
   \code{\link{dghyp}}
}
\examples{
QQplot(rnorm(1000),reference="normal"); 
QQplot(rexp(1000),reference="exp", rate=0.3); 
}
\keyword{hplot}

