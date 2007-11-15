\name{rACp}
\alias{rACp}
\title{ Simulate a Generalized Archimedean Copula representing p factors }
\description{
  Generate a sample from a generalized Archimedean copula.
}
\usage{
rACp(n, d, theta, A, name = "gumbel")
}
\arguments{
  \item{n}{ sample size of copula to be generated }
  \item{d}{ dimension of the simulated copula }
  \item{theta}{vector of parameters for each of p factors corresponding to underlying exchangeable copulas }
  \item{A}{ matrix with dimension d by p containing asymmetry parameters.  Rows must sum to 1. }
  \item{name}{ name of underlying exchangeable Archimedean copula from following list: 
     clayton, gumbel, frank, BB9, GIG. Names are case-sensitive.}
}
\details{
  See pages 224-6 of QRM for a bivariate example of this copula.  The idea carries naturally to higher dimensions.
}
\value{
  a matrix of dimension n times d where rows are realizations.
 
}
\seealso{ 
   \code{\link{rAC}},
   \code{\link{rcopula.Gumbel2Gp}}
   \code{\link{rcopula.GumbelNested}}
}
\examples{
alpha <- c(0.95,0.7);
wtmatrix <- cbind(alpha, 1-alpha);
data <- rACp(1000,d=2, theta=c(4,1), A=wtmatrix, name="gumbel");
}
\keyword{distribution }

