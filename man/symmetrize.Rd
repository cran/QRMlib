\name{symmetrize}
\alias{symmetrize}
\title{
  Ensure Symmetric Matrix 
}
\description{
ensures a matrix that should be symmetric is really symmetric 
}
\usage{
symmetrize(matrix)
}
\arguments{
\item{matrix}{
a matrix that should be symmetric 
}
}
\value{
a matrix that is symmetric 
}
\details{
deals with situations where rounding errors cause symmetric 
matrices to appear asymmetric 
}
\examples{
\dontrun{
#lines of code taken from fit.mst() in functionsNormix.R
# ...
Sigma <- var(data);
Sigma <- symmetrize(Sigma);
beta <- as.vector(solve(Sigma) %*% gamma);
mean <- as.numeric(mu+EW*gamma);
covariance <- EW*Sigma + varW*outer(gamma,gamma);
}
}
\keyword{array}

