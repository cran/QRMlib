\name{rBB9Mix}
\alias{rBB9Mix}
\title{
  Mixture Distribution Yielding BB9 Copula 
}
\description{
generates random sample from mixing distribution required for 
sampling from Joe's BB9 copula using Laplace transform method 
}
\usage{
rBB9Mix(n, theta)
}
\arguments{
\item{n}{
size of sample 
}
\item{theta}{
values of two parameters of BB9 copula; first must be positive; 
second must be greater than one. 
}
}
\value{
random sample of size n 
}
\details{
see page 224 of QRM. Algorithm essentially generates V 
corresponding to Joe's BB9 copula. For this copula 
algorithm uses fairly naive rejection and is SLOW! 
}
\section{References}{
Joe, H. Multivariate Models and Dependence Concepts, 
Chapman and Hall, 1997. See page 154 for BB9 copula. 
}
\seealso{
   \code{\link{rAC}}
}
\keyword{distribution}

