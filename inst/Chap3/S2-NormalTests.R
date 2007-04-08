# S-Plus script developed by Professor Alexander McNeil, mcneil@math.ethz.ch
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and DJ data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.2.1\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the DJ data set (the Dow Jones for 30 stocks) AND saved 
#the workspace, you may comment out the following line:
data(DJ);
#################################################


#####INITIAL CODE ################
#The initial code is repeated from S1-MomentEstimation.R in Chapter 3. Hence you
#may not need to run it again!!!!
### Some real multivariate return data.
# DJ is the Dow Jones 30 Components price timeSeries from 01/02/1991 - 12/29/2000. Convert 
# to a RETURN timeSeries.
Ret.DJ <- mk.returns(DJ);

#Choose only 10 of the 30 stocks:
selection1 <- c("AXP","EK","BA","C","KO","MSFT","HWP","INTC","JPM","DIS");
partialDJ30dailyTS <- Ret.DJ[,selection1];
#Choose only the data from 1/1/1993 to 12/31/2000.  Note 'from' date must be
#day prior to desired start date.
#In version 240.10068, fCalendar uses cut() rather than cutSeries() to select a subset from timeseries:
partialDJ30daily <- cut(partialDJ30dailyTS,from="1992-12-31", to="2000-12-31");
dim(partialDJ30daily@Data);
rm(partialDJ30dailyTS);
partialDJ30dailyMatrix <- seriesData(partialDJ30daily);
partialDJ30weeklyMatrix <- seriesData(aggregateWeeklySeries(partialDJ30daily, FUNC= colSums));
partialDJ30monthlyMatrix <- seriesData(aggregateMonthlySeries(partialDJ30daily, FUNC= colSums));
partialDJ30quarterlyMatrix <- seriesData(aggregateQuarterlySeries(partialDJ30daily, FUNC= colSums));

#########END INITIAL CODE##################

## Simulated data
ndata <- rmnorm(2000,rho=0.7,d=3); #from multivariate normal
tmp <- rcopula.clayton(2000,theta=1,d=5); #simulated from clayton copula
#apply the qnorm() function on the columns (=2 for 2nd parameter) of tmp.  This should produce
#the normal quantiles on each of the 5 dimensions of the copula
data.margnormal <- apply(tmp,2,qnorm);
#Create a matrix of scatterplots:
pairs(data.margnormal);
#generate random values from multivariate student-t distribution:
tdata <- rmt(2000,rho=0.7,d=10);

### Mardia's tests of multivariate normality for simulated data. 
#Called from functionsMultivariate.R.
MardiaTest(ndata);
MardiaTest(data.margnormal);
MardiaTest(tdata);

### Mardia's tests of multivariate normality: see pages 69-70 in QRM book
#Called from functionsMultivariate.R in QRMlib.
MardiaTest(partialDJ30dailyMatrix);
MardiaTest(partialDJ30weeklyMatrix);
MardiaTest(partialDJ30monthlyMatrix);
MardiaTest(partialDJ30quarterlyMatrix);

############## Univariate tests of normality (Jarque-Bera)
#No 'finmetrics' module in R: it exists only in S-Plus
#module(finmetrics)

#  ******IMPORTANT: must load a package ************
#The following R-language tests REQUIRE that fBasics package be loaded.  If the following
#library call fails, go to the CRAN site and download the fBasics package and the MASS
#package required by fBasics:
library(fBasics); #load fBasics package
#R's normal test does NOT have a parameter to remove NA from data:
#normalTest(data.margnormal,"jb", na.rm=T)
normalTest(data.margnormal,"jb");
normalTest(partialDJ30dailyMatrix, "jb");

#Two equivalent Jarque-Bera tests available in R include:
jbTest(data.margnormal); #high-precision
jarqueberaTest(data.margnormal);
weekly.pvals <- jbTest(partialDJ30weeklyMatrix);
weekly.pvals;
weekly.pvals <- normalTest(partialDJ30weeklyMatrix, "jb");
weekly.pvals;

#The normalTest() returns an fHTEST object with many slots (attributes) and corresponding
#values. Rather than returning all the values from the Jarque-Bera test, we can return only
# the p-values.  To see all the names of the many slots (or attributes)returned from 
#normalTest(),run unclass() on the returned object:
#  unclass(weekly.pvals)
# We then see that the p-value is returned from the "test" slot (attribute) so we can 
# follow the function invocation with @test$p.value to get just the p-value returned.
daily.pvals <- normalTest(partialDJ30dailyMatrix, "jb")@test$p.value;
daily.pvals;
weekly.pvals <- normalTest(partialDJ30weeklyMatrix, "jb")@test$p.value;
weekly.pvals;
quarterly.pvals <- normalTest(partialDJ30quarterlyMatrix, "jb")@test$p.value;
quarterly.pvals;
monthly.pvals <- normalTest(partialDJ30monthlyMatrix, "jb")@test$p.value;
monthly.pvals;
#Round off the pvalues to 3 decimal places:
base:::round(data.frame(daily.pvals,weekly.pvals,monthly.pvals,quarterly.pvals),3);


# Tests of normality based on Mahalanobis distance. Called from functionsMultivariate.R.
#Note the tests based on simulated data pass with reasonable K-S values:
jointnormalTest(ndata);
jointnormalTest(data.margnormal);
#Note the tests on the ten stocks selected from DJ30 fail the test miserably
#except possibly the quarterly values.  The QQ plots are very revealing. See p. 72 in QRM Book.
jointnormalTest(partialDJ30dailyMatrix);
jointnormalTest(partialDJ30weeklyMatrix);
jointnormalTest(partialDJ30quarterlyMatrix);
jointnormalTest(partialDJ30monthlyMatrix);

#Detach the library packages associated with fBasics that we loaded to run 
# normality tests (normalTest() and jarqueberaTest().
detach("package:MASS");
detach("package:fBasics");


