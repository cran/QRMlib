\name{rcopula.Gumbel2Gp}
\alias{rcopula.Gumbel2Gp}
\title{
  Gumbel Copula with Two-Group Structure 
}
\description{
generates sample from a Gumbel copula with two-group structure 
constructed using three Gumbel generators 
}
\usage{
rcopula.Gumbel2Gp(n=1000, gpsizes=c(2, 2), theta=c(2, 3, 5))
}
\arguments{
\item{n}{
sample size 
}
\item{gpsizes}{
vector of length two containing sizes of the groups 
}
\item{theta}{
parameter vector of length 3 giving parameters of the three Gumbel generators 
}
}
\value{
matrix of dimension n by sum(gpsizes) where rows are realizations 
}
\details{
see page 227 of QRM for an example of construction 
}
\seealso{
   \code{\link{rAC}},
   \code{\link{rcopula.AGumbel}},
   \code{\link{rcopula.GumbelNested}}
}
\examples{
data <- rcopula.Gumbel2Gp(n=3000,gpsizes=c(3,4),
                     theta=c(2,3,5)); 
pairs(data); 
}
\keyword{distribution}

