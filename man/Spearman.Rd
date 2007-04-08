\name{Spearman}
\alias{Spearman}
\title{
  Spearman's Rank Correlation 
}
\description{
calculates a matrix of Sperman's rank correlations 
}
\usage{
Spearman(data)
}
\arguments{
\item{data}{
data matrix 
}
}
\value{
matrix of rank correlations 
}
\details{
see pages 229-230 in QRM 
}
\seealso{
   \code{\link{Kendall}}
}
\examples{
data <- rmnorm(1000,d=3,rho=0.5); 
Spearman(data); 
}
\keyword{multivariate}

