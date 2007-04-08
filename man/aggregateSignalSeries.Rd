\name{aggregateSignalSeries}
\alias{aggregateSignalSeries}
\title{aggregateSignalSeries() method }
\description{
This is a substitute for the S-Plus language method\cr
\emph{aggregateSeries(signalSeries, FUN={max, mean, colAvg,...},\cr
       by=90,...)}. 
}
\usage{
aggregateSignalSeries(x, pos, AGGFUNC, together = FALSE, 
   drop.empty = TRUE, include.ends = FALSE, adj, offset, colnames, by)
}
\arguments{
  \item{x}{The data series to which the AGGFUNC will be applied}
  \item{pos}{ a numeric sequence describing the respective positions of each element in the data set}
  \item{AGGFUNC}{ the function to be applied to the data set x }
  \item{together}{if TRUE, pass all columns of x together into AGGFUNC; default is to pass each column 
separately into AGGFUNC for each aggregation block. 
 }
  \item{drop.empty}{logical value telling whether or not to drop aggregation blocks with no positions to 
aggregate from the output (default is to drop them) }
  \item{include.ends}{logical value telling whether to include positions before the first given 
aggregation block and after the last in the first/last blocks; default would not include those 
positions in the output at all. }
  \item{adj}{if provided, adjust the positions of the output series so that they lie a fraction adj 
towards the blocks ending position; default is to use the lower end of the block for the output 
position. 0.5 would use the center of the aggregation block for the output position, and 1.0 would 
use the upper end of the block. }
  \item{offset}{ as an alternative to adj, you can provide a constant offset to add to the lower 
end of the aggregation block to get the output series positions.}
  \item{colnames}{new column names for the output series. Default is to use the same names as the 
input series if the output series has the same width}
  \item{by}{The number of positions to include for each function application.  For example by=90 implies
the function will be applied to successive groups of 90 data items.}
}
\details{
  Input a signalSeries object as parameter x.  Input an a function (AGGFUNC) to apply to aggregate
data into many smaller subsamples.  Use either the \sQuote{pos} or \sQuote{by} parameter to indicated 
how to aggregate the data.  E.g. \sQuote{by=90} will chop the data into separate segments of length 90.  
The AGGFUNC will then be applied to each aggregation (segment).  The R-language 
\bold{aggregateSignalSeries()} function allows the use of a function evaluation (like Pearson or 
Kendall correlations) to create from data aggregated into granular group using the by parameter.  
E.g the  \emph{by=90} parameter will divide the dataset into groups of 90 observations and will 
apply the input function to each group of 90 data items.  Hence in 360 total observations, a total 
of four separate correlation functions may be evaluated on aggregated data sets each containing 
90 observations.
}
\value{
 A new signalSeries whose positions were adjusted via the \sQuote{by} parameter.  Hence the new 
signalSeries data slot contains types returned by the AGGFUN.  For example, if AGGFUNC is 
\sQuote{pearson} as in the example,then the data slot contains a vector of correlation coefficients  
each calculated by splitting the input data into successive blocks specified by using the number of 
items in \sQuote{by} for each new block.  For each block, the AGGFUNC is applied to each column (or  
all columns joined if parameter 'together=TRUE') to calculate the appropriate result.  The data slot 
contains the result applied to each successive block. 
}
\seealso{
   \code{\link{signalSeries}}
}
\examples{
\dontrun{
set.seed(13);
m <- 90;
n <- 3000;
#Generate a 'matrix' class of simulated values with 2 columns and m*n rows
dataSim <- rmt(m*n,df=3,rho=0.5,d=2);
#create a signal series from simulated data:
dataSimSS <- signalSeries(dataSim);
#local function
pearson <- function(x) cor(x)[1,2];
pearson.cors <- aggregateSignalSeries(dataSimSS,by=m,
                             together=T,AGGFUNC=pearson);
#Extract the data part only to see a vector of correlation 
#coefficients for each contiguous subblock
#in the entire original series.
pearson.cors.data <- pearson.cors@data;
}
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}

