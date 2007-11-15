# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# QRMlib 1.4.2
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and DJ data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the DJ timeSeries (the Dow Jones for 30 stocks) AND 
#saved the workspace, you may comment out the following line:
data(DJ);
#Alternatively, if you want to load the dataframe instead of timeSeries,
#activate the following line:
#data(DJ.df);
#################################################


ILLUSTRATION OF STANDARD METHODS.  See section 2.3 in QRM book.

# Get the risk factor Data (Z) in S-Plus
#DJ is a time series for the 30 industrials.  It is included in the data of the QRMLib.
#Through R-2.5.1, timeSeries class originally belong in package fCalendar. 
#Version 221.10065 used cutSeries()method to select data only between the 'to' and 'from' dates. 
#Version 240.10068 used cut(). Both required the day PRIOR to desired start in "from".
#Sdata <- cut(DJ, from="1992-12-30", to="2000-12-31");
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
Sdata <- window(DJ, from="1992-12-31", to="2000-12-31");
#Indicate the companies whose data you want in the reduced timeSeries:
tsSelections <- c("GE","INTC","KO","JNJ");
Sdata <- Sdata[,tsSelections];
#fCalendar 221.10065 used logSeries(). fCalendar 240.10068 used log(). fSeries 260.72 uses log()
Zdata <- log(Sdata); 

#In fCalendar 221.10065 and R-2.2.1, plot.timeSeries() worked to plot multiple series on same graph.
#plot(Zdata, ylab="log of timeSeries"); 
#Later versions no longer work with more than one column of data at a time.
#Hence we now call a new utility function from functionsUtility.R replacing plot.timeSeries():
#Plot all 4 columns on same graph (no need to use colvec parameter):
plotMultiTS(Zdata, reference.grid=TRUE);
#plot only columns 2 and 3 on the graph:
plotMultiTS(Zdata, colvec= c(2,3),reference.grid=TRUE, format="%Y-%m");

#Construct the risk factor returns (X)
Xdata <- mk.returns(Sdata);
#Convert from timeSeries to matrix:
X <- seriesData(Xdata);
#show a scatterplot of the four different returns for GE,INTC,KO,JNJ
pairs(X);


#McNeil extracts current prices from last row of Zdata series and then makes an exp() transform. 
#current.prices <- exp(Zdata[lastRow,]) #McNeil's original S-Plus way
#However, each element of Zdata is already a log() transform of Sdata.  Hence the current price may just as 
#well be removed as the last row of Sdata. Length of both Zdata and Sdata is length(Sdata@positions).  
#Extract last row (taking all columns) of Sdata time series to get 'current.prices' timeSeries as mentioned above. 
current.prices <- Sdata[length(Sdata@positions),];
current.prices; #display the current prices

#Let alpha represent the number of shares in each of the four assets of our portfolio; assume we have 1000 shares of each.
alpha <- c(1000,1000,1000,1000);
#Extract the current.prices as a vector from the timeSeries and transpose:
current.prices <- as.vector(t(seriesData(current.prices)));
#Calculate the portfolio value:
Value <- sum(current.prices * alpha);
#Calculate the corresponding portfolio weights:
portweights <- (current.prices * alpha)/Value;
portweights;

# Construct a 'loss operator function'.
loss.operator <- function(x,lop.weights,lop.value, linear=FALSE)
  # parameter x should be a matrix or vector of returns (risk-factors)
  # parameter lop.weights should be a vector of weights to apply to return vector
  # parameter lop.value should be a scalar representing a portfolio value
  # If risk-factor (e.g. returns) were calculated 
{
  if ((!(is.vector(x))) & (!(is.matrix(x)))) 
    stop("x must be vector or matrix with rows corresponding to risk factor return observations"); 
  nriskfactors <- length(lop.weights);
  ndatarows <- 1;
  if (is.matrix(x))
    ndatarows <- dim(x)[1];
  #Build a weight-matrix.  Let each row contain weights for the nriskfactors (e.g.four returns from 4 stocks)
  wtMatrix <- matrix(lop.weights,nrow=ndatarows,ncol=nriskfactors,byrow=TRUE);
  #Now multiply each element in the weight matrix by total portfolio value at start date. Hence each row of
  #tmpMatrix contains columns representing the starting individual risk factor values (e.g. four values for 4 stocks) 
  # Each row has currently identical values, the starting values for each individual risk factor
  tmp.matrix <- lop.value *wtMatrix;
  # tmp.matrix <- lop.value*matrix(lop.weights,nrow=ndatarows,ncol=nriskfactors,byrow=TRUE)
  #The default is nonlinear which means the function which calculated the return factor (e.g. mk.returns() used
  #a default log() (or some other nonlinear) transform.  Hence using exp(x) will reconvert the log() so the summand
  #will reflect actual changes in value  
  if (linear) 
     summand <- x*tmp.matrix
  else 
     summand <- (exp(x)-1)*tmp.matrix;
  #The following apply() function applies the 'sum' function to the ROWS of the summand matrix (since 2nd parameter is 1)
  #Hence we return a vector where each row is the change in value summed across each of the risk factors on a single day.
  #By taking the negative we convert to losses
  -as.vector(apply(summand,1,sum))
}

# Using this loss operator
loss <- loss.operator(X,portweights,Value); #losses on all 4 stocks in X throughout 2020 sample days
loss.operator(as.vector(X[1,]),portweights,Value); #losses on all 4 stocks in X on first day only
loss.operator(as.vector(X[1,]),portweights,Value,linear=TRUE); #losses on all 4 stocks in X on first day only--linearized



#Implement variance-covariance analysis
muhat <- as.vector(apply(X,2,mean));
#R-language does NOT have an unbiased=FALSE option.  Hence we must multiply result by (n-1)/n
#Sigmahat  <- var(X,unbiased=F) #S-Plus biased version
number <- dim(X)[1];
Sigmahat  <- var(X)*(number-1)/number;  #R biased version

meanloss <- -sum(portweights*muhat)*Value
varloss <- Value^2 *(t(portweights) %*% Sigmahat %*% portweights);
VaR99 <- meanloss + sqrt(varloss)*qnorm(0.99);
ES99 <- meanloss +sqrt(varloss)*dnorm(qnorm(0.99))/0.01;

#Implement a historical simulation analysis 
hsdata <- loss.operator(X,portweights,Value);
hist(hsdata,nclass=20);
qqnorm(hsdata);
VaR99.hs <- quantile(hsdata,0.99);
ES99.hs <- mean(hsdata[hsdata > VaR99.hs]);
#R-language requires seriesPositions() rather than positions()
hsdata.ts <- timeSeries(hsdata,seriesPositions(Xdata));
#Beginning with fCalendar 240.10068,  you must be sure to pass the type to the plot():
plot(hsdata.ts,type="l");

#Implement a Monte Carlo simulation analysis
#You might import simulated data from S-Plus to test results:
X.new <- rmnorm(10000,Sigma=Sigmahat,mu=muhat);
mcdata <- loss.operator(X.new,portweights,Value);
hist(mcdata);
qqnorm(mcdata);
VaR99.mc <- quantile(mcdata,0.99);
ES99.mc <- mean(mcdata[mcdata > VaR99.mc]);

#Draw pictures
hist(hsdata,nclass=20,prob=TRUE, xlab="HS data");
title(main="HS Loss Distribution");
abline(v=c(VaR99,ES99));
abline(v=c(VaR99.hs,ES99.hs),col=2,lty=2);
abline(v=c(VaR99.mc,ES99.mc),col=3,lty=3);

#Implement alternative Monte Carlo simulation analysis based on t
#Fit multivariate t to X data to get parameters for t; then simulate from t using those parameters:
model <- fit.mst(X);
X.newt <- rmt(10000,df=model$nu,Sigma=model$Sigma,mu=model$mu);
mcdatat <- loss.operator(X.newt,portweights,Value);
VaR99.mct <- quantile(mcdatat,0.99);
ES99.mct <- mean(mcdatat[mcdatat > VaR99.mct]);

# Add details to graph
abline(v=c(VaR99.mct,ES99.mct),col=4, lty=4);
box();
legend(-15000,0.00015, c("Variance Covariance","Historical Simulation","Monte Carlo (Gaussian)","Monte Carlo (t)"),col=c(1,2,3,4),lty=c(1,2,3,4));
