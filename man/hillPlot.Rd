\name{hillPlot}
\alias{hillPlot}
\title{ Create Hill Plot }
\description{
  Plot the Hill estimate of the tail index of heavy-tailed data, 
  or of an associated quantile estimate.
}
\usage{
hillPlot(data, option = c("alpha", "xi", "quantile"), start = 15, 
   end = NA, reverse = FALSE, 
   p = NA, ci = 0.95, auto.scale = TRUE, labels = TRUE, ...)
}
\arguments{
  \item{data}{ data vector}
  \item{option}{whether "alpha", "xi" (1/alpha) or "quantile" (a quantile estimate) should be plotted}
  \item{start}{lowest number of order statistics at which to plot a point}
  \item{end}{highest number of order statistics at which to plot a point}
  \item{reverse}{whether plot is to be by increasing threshold (TRUE) or increasing number of order statistics (FALSE)}
  \item{p}{probability required when option "quantile" is chosen}
  \item{ci}{probability for asymptotic confidence band; for no confidence band set ci to zero}
  \item{auto.scale}{whether or not plot should be automatically scaled; if not, xlim and ylim 
graphical parameters may be entered}
  \item{labels}{whether or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\details{
  This plot is usually calculated from the alpha perspective. For a generalized Pareto analysis of heavy-tailed 
  data using the gpd function, it helps to plot the Hill estimates for xi. 
  See pp. 286-289 in QRM.  Especially note that Example 7.28 suggests the best estimates
  occur when the threshold is very small, perhaps 0.1%-0.5% of the sample size (10-50 order
  statistics in a sample of size 1000.  Hence you should NOT be using a 95% threshold for Hill 
  estimates.
}
\value{
  None
}

\seealso{
   \code{\link{xiplot}},
   \code{\link{plotTail}}
}

\examples{
data(danish);
#Run hillPlot to show what happens with the Hill Plot.  
#See Example 7.27, p. 287 in QRM
hillPlot(danish, option = "alpha", start = 5, end = 250, p = 0.99);
hillPlot(danish, option = "alpha", start = 5, end = 60, p = 0.99);
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{hplot}

