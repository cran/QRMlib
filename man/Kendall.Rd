\name{Kendall}
\alias{Kendall}
\title{
  Kendall's Rank Correlation 
}
\description{
calculates a matrix of Kendall's rank correlations 
}
\usage{
Kendall(data, noforce=TRUE)
}
\arguments{
\item{data}{
data matrix 
}
\item{noforce}{
must be set to FALSE if you really want to calculate Kendall's rank correlations 
for more than 5000 data (which will be slow) 
}
}
\value{
matrix of rank correlations 
}
\details{
see pages 229-230 in QRM 
}
\seealso{
   \code{\link{Spearman}},
   \code{\link{fit.tcopula.rank}}
}

\examples{
data <- rmnorm(1000,d=3,rho=0.5); 
Kendall(data); 
}
\keyword{multivariate}
\keyword{methods}

