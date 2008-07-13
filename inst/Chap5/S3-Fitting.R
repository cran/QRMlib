# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# QRMlib 1.4.4
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and ftse100 data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the ftse100, smi,and FXGBP.RAW timeSeries AND saved 
#the workspace, you may comment out the following lines:
data(ftse100);
data(smi);
data(FXGBP.RAW);
#Alternatively, if you want to load the dataframes instead of timeSeries,
#activate the following lines:
#data(ftse100.df);
#data(smi.df);
#data(FXGBP.RAW.df);
#################################################

#### COPULA FITTING
# Fitting copulas to data the R-way
#requires(fSeries) #version 260.72

#In R (unlike S-Plus), two time series being merged must have the same number of rows.
#The smi timeSeries starts on 11/09/1990 and the ftse100 starts on 01/02/1990. Both end 03/25/2004.
#Hence create a partial ftse100 time series by cutting all values between 1/2/1990-11/8/1990. so we 
#will use the remaining (cut) data from 1990-11-08 to 2004-03-25. Note the first date 1990-11-08 is
#excluded
#Through R-2.5.1, timeSeries class originally belong in package fCalendar. 
#Version 221.10065 used cutSeries()method to select data only between the 'to' and 'from' dates. 
#Version 240.10068 used cut(). Both required the day PRIOR to desired start in "from".
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
TS1 <- window(ftse100, "1990-11-09", "2004-03-25"); # TS1 is ftse100 truncated.

#If we try to merge the two series which have the same dates, we would hope to get two series of
#equal length. However, if we try to merge the series (use only the @Data slot from the 2nd series)
#we get an error message that series are NOT of same length. 
#INDEXES.RAW <- mergeSeries(TS3, smi@Data)
# Investigation shows that TS3 has 3379 entries while TS2 has 3331 entries.  TS2 is missing, 
#for example, dates 1990-12-31 and 1991-01-02.
#dim(TS1@Data)[1]; #shows 3379
#dim(smi@Data)[1]; #shows 3331
#When S-Plus finds data mismatches, its ts.union() function uses the 'before' method to copy data
#from the prior date to any missing date.  This is not done in R.  Hence we have to run a 2nd
#set of functions to align the data using all days available.  If a date has missing data, the
#"before" argument says to copy the data from the day before to the date with missing data.
TS1Augment <- alignDailySeries(TS1, method="before"); #gives 3490 observations
TS2Augment <- alignDailySeries(smi, method="before"); #gives 3490 observations
#Merge the two time series with ftse100 first. In fCalendar versions prior to 240.10068, we used
#mergeSeries().  In fCalendar 240.10068, mergeSeries() has been deprecated and we must use merge();
#INDEXES.RAW <- mergeSeries(TS1Augment, TS2Augment@Data);
INDEXES.RAW <- merge(TS1Augment,TS2Augment);
#clean up:
rm(TS1, TS1Augment, TS2Augment);
#Prior to R-2.5.0 we called plot.timeSeriesIts()from utilityFunction.R. However in R-2.5.0 both 
#'its' and 'base' have a method called "names" and  R-2.5.0 cannot distinguish them properly. 
#"names" when both 'base' and 'its' are always loaded. 
#Hence we now call a new utility function from functionsUtility.R replacing plot.timeSeriesIts():
#Plot both columns on same graph (no need to use colvec parameter):
plotMultiTS(INDEXES.RAW, reference.grid=TRUE);

#Create a time series of returns from the prices.  mk.returns is in functionsUtility.R.
#By not passing a type="relative" argument, we are implicitly using a type="log" argument
# and generating logarithmic returns.
INDEXES <- mk.returns(INDEXES.RAW);
#As above call new utilityFunction.R plotMultiTS():
plotMultiTS(INDEXES);


#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
PARTIALINDEXES <- window(INDEXES, "1994-01-01", "2003-12-31"); 

#Now create a data matrix from the just-created timeSeries 
data <- seriesData(PARTIALINDEXES);
#Keep only the data items which are non-zero for both smi and ftse100
data <- data[data[,1]!=0 & data[,2] !=0,];
plot(data);

# Construct pseudo copula data. The 2nd parameter is MARGIN=2 
#when applying to columns and 1 applied to rows. Hence this says to
#apply the 'edf()' empirical distribtion function() to the columns
#of the data. 
Udata <- apply(data,2,edf,adjust=1)
plot(Udata)

# Compare various bivariate models. Fit gauss copula first
mod.gauss <- fit.gausscopula(Udata)
mod.gauss

#Calculate Spearman rank correlation for pseudo-copula data:
Spearman(Udata)

#Fit 2-dimensional Archimedian copula: choices are gumbel or clayton
#using pseudo data generated via edf() from observed data:
mod.gumbel <- fit.Archcopula2d(Udata,"gumbel")
mod.clayton <- fit.Archcopula2d(Udata,"clayton")

#Fit a t-copula to the pseudo-copula data generated via edf() from real data: 
mod.t <- fit.tcopula(Udata)
mod.t

sin(pi*Kendall(Udata)/2)
c(mod.gauss$ll.max, mod.t$ll.max, mod.gumbel$ll.max, mod.clayton$ll.max)


# Multivariate Fitting with Gauss and t: Simulation
# Create an equicorrelation matrix:
P <- equicorr(3,0.6)

#Simulate data for a Gaussian copula with Sigma=P
set.seed(117);
#The data will have 1000 rows and 3 columns since P is 3x3:
Udatasim <- rcopula.gauss(1000,Sigma=P);
#ALTERNTATIVELY, copy data simulated from S-Plus for usage
#Try to fit a copula to the simulated data. It should match.
mod.gauss <- fit.gausscopula(Udatasim);
mod.gauss;

#Calculate a matrix of Spearman Rank Correlations
Pstar <- Spearman(Udatasim);
#Evalutate the density of the Gauss copula using dcopula.gauss() and the matrix
# of Spearman Rank Correlations. This returns a vector. 
# Sum the vector elements.  See pp. 197 and 234 in QRM
sum(dcopula.gauss(Udatasim,Pstar,logvalue=TRUE));

#Reset the seed to a new value:
set.seed(113);
#Generate a new set of random data from a t-copula (10df) with the same Sigma matrix:
Udatasim2 <- rcopula.t(1000,df=10,Sigma=P)
#Fit the copula to the simulated data
mod.t <- fit.tcopula(Udatasim2);
mod.t;

#Now refit the copula to the same simulated data using (Kendall) rank correlations
#and the fit.tcopula.rank() method:
mod.t2 <- fit.tcopula.rank(Udatasim2);
mod.t2;
#Now fit the gauss copula to the same simulated data
mod.gauss <- fit.gausscopula(Udatasim2);
mod.gauss;

# ******Now use real FX data with Great Britain's Pound. 
#Once again, use the new utility function for R-2.5.0:
plotMultiTS(FXGBP.RAW, reference.grid=TRUE);
#Create a return time series; it will be logarithmic unless we pass argument type="relative"
tsretFXGBP <- mk.returns(FXGBP.RAW);
#In version 240.10068, fCalendar uses cut() rather than cutSeries() to select a subset from timeseries:
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
tsretFXGBP <- window(tsretFXGBP,from= "01/01/1997", to= "12/31/2003");
#Extract just the data, removing the date; effectively creates a data matrix
data <- seriesData(tsretFXGBP);

#Use the sum() function across rows (indicated by MARGIN=1)of the matrix.  If any row
#has missing values, an NA will be returned.  Hence this reports on which rows have
#missing values since there are 4 different exchange rates per row.  Missing data
#will be a row value of TRUE in 'missing'; otherwise row values will be FALSE.
missing <- is.na(apply(data,1,sum));
#REMOVE any rows with missing values, i.e. keep only data which is not missing values
data <- data[!missing,];
#print scatter plots of all data pairs (e.g GBP/USD,  GBP/EUR, GBP/JPY, GBP/CHF)
#remember-this is the real data
pairs(data);

# Construct pseudo copula data
#Apply the edf() function to the columns.  This will calculate the empirical cumulative
#probability for each observation in the time series.
#The edf() is the empirical distribution function. edf() takes an argument 'adjust' 
#which says whether the denominator should be 'n+1'(adjust=1) or not (adjust=0). The 
#default is 0 (F). Hence the resulting Udata will represent the value of the ecdf 
#(empirical cdf) represented by each return.
#The 2nd argument of apply is the MARGIN: 1 = rows, 2=columns, c(1,2) for both rows 
#and columns)
UdataFX <- apply(data,2,edf,adjust=TRUE);
#print scatter plots of all data pairs (e.g GBP/USD,  GBP/EUR, GBP/JPY, GBP/CHF) 
#generated from Empirical Distribution Function edf
pairs(UdataFX);
# Fit copulas using rank correlations to the ECDF data.
mod.t <- fit.tcopula.rank(UdataFX);
mod.t;
#Get Spearman rank correlations
Pstar <- Spearman(UdataFX);
Pstar;
sum(dcopula.gauss(UdataFX,Pstar,logvalue=TRUE));
Pstar;
