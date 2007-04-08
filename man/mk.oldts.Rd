\name{mk.oldts}
\alias{mk.oldts}
\title{
  Make Old Style Time Series 
}
\description{
This function converts a 'timeSeries' object (from package fCalendar in R) 
to an 'its' object (an irregular times series from package 'its' in R).  
The conversion is useful when plotting irregular (not-evenly spaced in time) 
series. 
}
\usage{
mk.oldts(timeseries)
}
\arguments{
\item{timeseries}{
a timeSeries object 
}
}
\value{
an old-style 'its' object 
}
\details{
this is occasionally useful for plotting 
}
\keyword{utilities}

