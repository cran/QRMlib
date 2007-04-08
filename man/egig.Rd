\name{EGIG}
\alias{EGIG}
\title{Estimate Moments of GIG Distribution}
\description{
Calculates moments of univariate generalized inverse Gaussian (GIG) distribution 
}
\details{
Normal variance mixtures are frequently obtained by perturbing the variance component
of a normal distribution; we multiply by the square root of a mixing variable assumed 
to have a GIG (generalized inverse gaussian) distribution depending upon three parameters
lambda, chi, and psi. See p.77 in QRM.\cr

Normal mean-variance mixtures are created from normal variance mixtures by applying another
perturbation of the same mixing variable to the mean component of a normal distribution.
These perturbations create Generalized Hyperbolic Distributions.  See pp. 78-81 in QRM.

Also see page 497 of QRM Book for a description of the GIG distribution.

}
\usage{EGIG(lambda, chi, psi, k=1)}
\arguments{
\item{lambda}{lambda parameter}
\item{chi}{chi parameter}
\item{psi}{psi parameter}
\item{k}{order of moment}
}
\value{mean of distribution}
\seealso{
   \code{\link{rGIG}}
   \code{\link{ElogGIG}}
}

\keyword{distribution}

