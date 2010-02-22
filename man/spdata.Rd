\name{spdata}
\docType{data}
\alias{spdata}
\title{
  Standard and Poors Default Data 
}
\description{
The \code{spdata} timeSeries dataset has 100 rows and 3 columns. 
It contains default data for A, BBB, BB, B and C-rated companies for the years 1981 to 2000 
}
\usage{
data(spdata)
}
\format{
a matrix containing 100 rows and 4 columns.\cr
The colums are: \cr 
\tabular{ll}{
\emph{rating} \tab rating category (A, BBB, BB, B, CCC) \cr
\emph{firms} \tab number of companies in rating category \cr
\emph{number of defaults} \tab number of companies defaulting in category}
The rows are the years from 1981-2000
}

\source{
Standard and Poors Credit Monitor 
}
\seealso{
   \code{\link{spdata.df}},
   \code{\link{spdata.raw}},
   \code{\link{momest}}
}

\examples{
#Must attach MASS and nlme libraries to run mixed effect regression model.
#library(MASS) is attached by fSeries which is a dependency of QRMlib
#and hence is already loaded
library(nlme);
#Load timeSeries:
data(spdata); #timeSeries
ratingval <- spdata@recordIDs$rating;
yearval <- as.numeric(spdata@recordIDs$DATE);
#Use R- library MASS to get glmmPQL which runs a mixed-effects model. 
#It will measure random effects and fixed effects.
#'year' -'ratings' determine the unique results(20 years 1981-2000 with 5 obligor 
#class ratings each year)
results <- glmmPQL(cbind(defaults,firms-defaults) ~ -1 + ratingval, 
     random = ~1| yearval, family=binomial(probit), data=spdata);
results;
summary(results);
summary(results)$tTable[,1];
detach("package:nlme");
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{datasets}

