\name{fit.AC}
\alias{fit.AC}
\title{ Fit Archimedean Copula }
\description{
  Fit an Archimedean copula via maximum likelihood
}
\usage{
fit.AC(Udata, name = "clayton")
}
\arguments{
  \item{Udata}{ matrix of copula data consisting of n rows of d-dimensional vector observations}
  \item{name}{ name of copula, e.g. "clayton" or "gumbel" }
}
\details{
  see documentation of dcopula.AC for information on extending.
This function can be used in place of the older function fit.Archcopula2d().
}
\value{
  list containing parameter estimate, standard error, value of log-likelihood at maximum and convergence flag
}

\seealso{\link{dcopula.AC}}
\keyword{ distribution }

