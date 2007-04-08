\name{Pdeconstruct}
\alias{Pdeconstruct}
\title{
  Disassemble a Correlation Matrix for ML Copula Fitting 
}
\description{
takes a correlation matrix P and returns the elements of a 
lower-triangular matrix A with ones on the diagonal such that 
P is the corelation matrix corresponding to the covariance matrix AA' 
}
\usage{
Pdeconstruct(P)
}
\arguments{
\item{P}{
a correlation matrix 
}
}
\value{
elements of a lower-triangular matrix 
}
\details{
see page 235 in QRM 
}
\seealso{
   \code{\link{Pconstruct}},
   \code{\link{fit.gausscopula}},
   \code{\link{fit.tcopula}}
}
\examples{
P <- Pconstruct(c(1,2,3,4,5,6)); 
Pdeconstruct(P); 
}
\keyword{array}
