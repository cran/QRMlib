\name{jointnormalTest}
\alias{jointnormalTest}
\title{
  Test of Multivariate Normality 
}
\description{
provides test of multivariate normality based on analysing Mahalanobis distances 
}
\usage{
jointnormalTest(data, dist="chisquare")
}
\arguments{
\item{data}{
matrix of data with each row representing an observation 
}
\item{dist}{
"chisquare" performs test against chi-squared distribution, 
which is an approximation; "beta" performs test against a scaled beta 
}
}
\value{
p-value for Kolmogorov-Smirnov test 
}
\section{Side Effects}{
a QQplot against the reference distribution is created 
}
\details{
see pages 69-70 of QRM 
}
\seealso{
   \code{\link{MardiaTest}}
}
\examples{
data(DJ);
Ret.DJ <- mk.returns(DJ);
selection1 <- c("AXP","EK","BA","C","KO","MSFT",
                "HWP","INTC","JPM","DIS");
partialDJ30dailyTS <- Ret.DJ[,selection1];
#Choose only the data from 1/1/1993 to 12/31/2000.  
partialDJ30daily <- window(partialDJ30dailyTS,from="1993-01-01", 
                               to="2000-12-31");
partialDJ30dailyMatrix <- seriesData(partialDJ30daily);
#Note the tests on the ten stocks selected from DJ30 fail the test miserably
#except possibly the quarterly values.  The QQ plots are very revealing. 
#See p. 72 in QRM Book.
jointnormalTest(partialDJ30dailyMatrix);
}
\keyword{multivariate}

