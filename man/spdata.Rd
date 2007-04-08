\name{spdata}
\docType{data}
\alias{spdata}
\title{
  Standard and Poors Default Data 
}
\description{
The \code{spdata} data frame has 100 rows and 4 columns. 
It contains default data for A, BBB, BB, B and C-rated companies for the years 1981 to 2000 
}
\usage{
data(spdata)
}
\format{
a matrix containing 100 rows and 4 columns.\cr
The colums are: \cr 
\tabular{ll}{
\emph{year} \tab year of default \cr
\emph{rating} \tab rating category (A, BBB, BB, B, CCC) \cr
\emph{firms} \tab number of companies in rating category \cr
\emph{number of defaults} \tab number of companies defaulting in category}

The rows are the years from 1981-2000
}

\source{
Standard and Poors Credit Monitor 
}
\seealso{
   \code{\link{spdata.raw}},
   \code{\link{momest}}
}

\examples{
\dontrun{
library(MASS);
library(nlme);
#the data being used is spdata which has 100 rows and 4 columns: 
#'year', 'rating', 'firms', defaults'
#Use R- library MASS to get glmmPQL
#'year' -'ratings' determine the unique results(20 years 1981-2000 with 5 obligor 
#class ratings each year
results <- glmmPQL(cbind(defaults,firms-defaults) ~ -1 + rating, 
    random = ~1| year, family=binomial(probit), data=spdata);
summary(results);
summary(results)$tTable[,1];
detach("package:nlme");
detach("package:MASS");
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{datasets}

