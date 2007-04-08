\name{rcopula.AGumbel}
\alias{rcopula.AGumbel}
\title{
  Generate Asymmetric Gumbel Copula 
}
\description{
generates sample from asymmetric Gumbel copula 
}
\usage{
rcopula.AGumbel(n, theta, alpha=rep(1, d), d=2)
}
\arguments{
\item{n}{
sample size 
}
\item{theta}{
scalar parameter of an exchangeable Gumbel copula 
}
\item{alpha}{
vector of length d containing parameters of asymmetry 
}
\item{d}{
dimension of copula 
}
}
\value{
a matrix of dimension n times d where rows are realizations 
}
\details{
see pages 224-226 of QRM for bivariate example of this copula; 
the idea obviously carries over to higher dimensions 
}
\seealso{
   \code{\link{rAC}},
   \code{\link{rcopula.Gumbel2Gp}},
   \code{\link{rcopula.GumbelNested}}
}
\examples{
data.AGumbel <- rcopula.AGumbel(10000,theta=4,
                     alpha=c(0.95,0.7)); 
plot(data.AGumbel); 
data.AGumbel <- rcopula.AGumbel(5000,theta=2,
                  c(0.1,0.7,0.8,0.9)); 
pairs(data.AGumbel); 
hist(data.AGumbel[,4]);
}
\keyword{distribution}

