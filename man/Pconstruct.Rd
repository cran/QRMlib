\name{Pconstruct}
\alias{Pconstruct}
\title{
  Assemble a Correlation Matrix for ML Copula Fitting 
}
\description{
takes a vector of values representing the terms of a lower 
triangular matrix A with ones on the diagonal and calculates 
the correlation matrix corresponding to the covariance matrix AA' 
}
\usage{
Pconstruct(theta)
}
\arguments{
\item{theta}{
elements of a lower triangular matrix A with ones on the diagonal 
}
}
\value{
correlation matrix corresponding to cobariance matrix AA' 
}
\details{
see page 235 in QRM 
}
\seealso{
   \code{\link{Pdeconstruct}},
   \code{\link{fit.gausscopula}},
   \code{\link{fit.tcopula}}
}
\examples{
P <- Pconstruct(c(1,2,3,4,5,6)); 
\dontrun{
eigen(P); 
Pdeconstruct(P); 
}
}
\keyword{array}

