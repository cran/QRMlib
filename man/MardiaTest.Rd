\name{MardiaTest}
\alias{MardiaTest}
\title{
  Mardia's Tests of Multinormality 
}
\description{
conducts Mardia's tests of multinormality based on multivariate 
skewness and kurtosis statistics 
}
\usage{
MardiaTest(data)
}
\arguments{
\item{data}{
data matrix 
}
}
\value{
vector consisting of skewness statistic, p-value for skewness statistic, 
kurtosis statistic and p-value for kurtosis statistic 
}
\details{
see page 70 of QRM 
}
\seealso{
   \code{\link{jointnormalTest}}
}

\examples{
data <- rmnorm(1000,d=10,rho=0.6); 
MardiaTest(data); 
}
\keyword{multivariate}

