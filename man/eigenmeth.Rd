\name{eigenmeth}
\alias{eigenmeth}
\title{
  Make Matrix Positive Definite 
}
\description{
adjusts a negative definite symmetric matrix to make it positive definite 
}
\usage{
eigenmeth(mat, delta=0.001)
}
\arguments{
\item{mat}{
a symmetric matrix 
}
\item{delta}{
new size of smallest eigenvalues 
}
}
\value{
a positive-definite matrix 
}
\details{
see page 231 of QRM 
}
\seealso{
   \code{\link{fit.tcopula.rank}}
}

\keyword{array}

