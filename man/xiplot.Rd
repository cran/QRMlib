\name{xiplot}
\alias{xiplot}
\title{
  GPD Shape Parameter Plot 
}
\description{
creates a plot showing how the estimate of shape varies 
with threshold or number of extremes. 
}
\usage{
xiplot(data, models=30., start=15., end=500., reverse=TRUE, 
ci=0.95, auto.scale=TRUE, labels=TRUE, table=FALSE, ...)
}
\arguments{
\item{data}{
vector or time series of data 
}
\item{models}{
number of consecutive gpd models to be fitted; i.e. the number of different
thresholds at which to re-estimate xi; this many xi estimates will be plotted 
}
\item{start}{
lowest number of exceedances to be considered 
}
\item{end}{
maximum number of exceedances to be considered 
}
\item{reverse}{
should plot be by increasing threshold (TRUE) or number of extremes (FALSE) 
}
\item{ci}{
probability for asymptotic confidence band; for no confidence 
band set to FALSE 
}
\item{auto.scale}{
whether or not plot should be automatically scaled; if not, 
xlim and ylim graphical parameters may be entered 
}
\item{labels}{
whether or not axes should be labelled; default is TRUE 
}
\item{table}{
should a table of results be printed; default is FALSE 
}
\item{\dots}{
further parameters of xiplot function 
}
}
\details{
For every model "fit.GPD" is called. Evaluation may be slow. 
}
\seealso{
   \code{\link{fit.GPD}},
   \code{\link{MEplot}}
}
\examples{
# Shape plot of heavy-tailed Danish fire insurance data:
data(danish);
xiplot(danish); 
}
\keyword{hplot}
