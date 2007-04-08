# S-Plus script developed by Professor Alexander McNeil, mcneil@math.ethz.ch
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and nasdaq index daily data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.2.1\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the daily nasdaq index data set  AND saved 
#the workspace, you may comment out the following line:
data(nasdaq);
#################################################


# Analysis of Block Maxima with GEV
#These use timeSeries class requiring
#require(fCalendar)
#be loaded.

plot(nasdaq, type="l");
#The following plot includes gridlines:
#plot.timeSeriesIts(nasdaq);

nreturns <- -mk.returns(nasdaq);
plot(nreturns, type="l",ylab="negative returns");
#plot.timeSeriesIts(nreturns); #includes gridlines

#Original S-Plus code uses call to S-Plus function:
#     monthly.maxima <- aggregateSeries(nreturns,FUN=max,by="months")
#This function is actually a call into a C++/C DLL (S.dll).  A function named
#aggregateSeries() supposedly belongs to timeSeries class present in fSeries 
#for R but DOES NOT EXIST.  Hence two new R-language functions have been added to the
#functionsUtility.R code to aggregate daily data series into monthly and 
#quarterly data for R-language.
#Use the following customized functions from functionsUtility.R
monthly.maxima <- aggregateMonthlySeries(nreturns,FUN=max);
quarterly.maxima <- aggregateQuarterlySeries(nreturns, FUNC=max);
plot(monthly.maxima, type="l", main="NASDAQ Monthly Maxima",ylab="rates of return");
plot(quarterly.maxima, type="l", main="NASDAQ Quarterly Maxima",ylab="rates of return");

monthly.maxima <- seriesData(monthly.maxima);
mod1 <- fit.GEV(monthly.maxima);
quarterly.maxima <- seriesData(quarterly.maxima);
mod2 <- fit.GEV(quarterly.maxima);

# Estimate 40 quarter return level
qGEV(1-1/40,xi=mod2$par.ests[1],mu=mod2$par.ests[3],sigma=mod2$par.ests[2]);

# Estimate probability of new record in a quarter
1-pGEV(max(quarterly.maxima),xi=mod2$par.ests[1],mu=mod2$par.ests[3],sigma=mod2$par.ests[2]);

