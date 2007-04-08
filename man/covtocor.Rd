\name{CovToCor}
\alias{CovToCor}
\title{
  Covariance To Correlation Matrix 
}
\description{
extracts the correlation matrix from a covariance matrix 
}
\details{
This is a custom function built by Alexander McNeil.  It provides the same
functionality as R's built in cov2cor() method in the stats package.
}
\usage{
CovToCor(mat)
}
\arguments{
\item{mat}{
a covariance matrix 
}
}
\value{
a correlation matrix 
}
\seealso{
   \code{\link{equicorr}}
   \code{\link{cov2cor}}
}

\keyword{array}

