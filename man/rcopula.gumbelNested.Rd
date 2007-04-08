\name{rcopula.GumbelNested}
\alias{rcopula.GumbelNested}
\title{
  Gumbel Copula with Nested Structure 
}
\description{
generates sample from a d dimensional Gumbel copula with 
nested structure constructed using (d-1) Gumbel generators 
}
\usage{
rcopula.GumbelNested(n, theta)
}
\arguments{
\item{n}{
sample size 
}
\item{theta}{
vector of admissable Gumbel copula parameters of length (
d-1) ordered by increasing size 
}
}
\value{
matrix of dimension n by d where rows are realizations 
}
\details{
see page 226 of QRM for trivial thrivariate example 
}
\seealso{
   \code{\link{rAC}},
   \code{\link{rcopula.AGumbel}},
   \code{\link{rcopula.Gumbel2Gp}}
}
\examples{
data <- rcopula.GumbelNested(n=3000,theta=1:6); 
pairs(data); 
}
\keyword{distribution}

