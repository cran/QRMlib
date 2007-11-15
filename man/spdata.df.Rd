\name{spdata.df}
\docType{data}
\alias{spdata.df}
\title{
  Standard and Poors Default Data 
}
\description{
The \code{spdata.df} data.frame has 100 rows and 4 columns. 
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
   \code{\link{spdata}},
   \code{\link{spdata.raw}},
   \code{\link{momest}}
}

\examples{
#Easier to use data.frame than timeSeries for this set
data(spdata.df);
#the data being used is spdata.df which has 100 rows and 4 columns: 
#'year', 'rating', 'firms', defaults'
#Must attach MASS and nlme libraries to run mixed effect regression model
library(MASS);
library(nlme);
#Use R- library MASS to get glmmPQL which runs a mixed-effects model. 
#It will measure random effects and fixed effects.
#'year' -'ratings' determine the unique results(20 years 1981-2000 with 5 obligor 
#class ratings each year)
results <- glmmPQL(cbind(defaults,firms-defaults) ~ -1 + rating, 
    random = ~1| year, family=binomial(probit), data=spdata.df);
results;
summary(results);
summary(results)$tTable[,1];
rm(spdata.df);
detach("package:nlme");
detach("package:MASS");
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{datasets}

