\name{rAC}
\alias{rAC}
\title{
  Generate Archimedean Copula 
}
\description{
generates data from a multivariate Archimedean copula with
 arbitrary dimension using the mixture construction of 
Marshall and Olkin 
}
\usage{
rAC(name, n, d, theta)
}
\arguments{
\item{name}{
name of the Archimedean copula from following list: 
clayton, gumbel, frank, BB9 
}
\item{n}{
number of realizations 
}
\item{d}{
dimension of copula 
}
\item{theta}{
parameter(s) of copula 
}
}
\value{
a matrix of dimension n times d where rows are realizations 
}
\details{
this function may easily be augmented with further Archimedean copulas 
}
\seealso{
   \code{\link{rcopula.clayton}},
   \code{\link{rcopula.frank}},
   \code{\link{rcopula.gumbel}},
   \code{\link{rcopula.gauss}},
   \code{\link{rcopula.t}}
}

\examples{
#simulate values from Archimedan copula of type gumbel
rAC("gumbel", n=3000, d=7,theta =3); 
}
\keyword{distribution}

