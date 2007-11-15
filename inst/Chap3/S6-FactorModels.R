# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)

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

#if you have previously opened the dji data set (the dow jones index) AND saved the
#workspace, you may comment out the following line:
data(dji);
#Alternatively, if you want to load the dataframe instead of timeSeries,
#activate the following line:
#data(dji.df);
#################################################


#######################################################
#DIMENSION REDUCTION AND FACTOR MODELS
######################################################



######## Factor Models#######################

# Construct appropriate return data from DJ timeSeries object
#This time use relative returns rather than the default log returns:
tsRelPctRet.DJ <- 100*mk.returns(DJ,type="relative");
selection <- c("MO","KO","EK","HWP","INTC","MSFT","IBM","MCD","WMT","DIS")
tsRelPctRet.DJ <- tsRelPctRet.DJ[,selection];
#In version 240.10068, fCalendar uses cut() rather than cutSeries() to select a subset from timeseries:
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
Xdata.ts <- window(tsRelPctRet.DJ, from="1992-01-01", to="1998-12-30");
#Extract the data only:
Xdata <- seriesData(Xdata.ts);
#produce a matrix of scatterplots:
pairs(Xdata);

tsRelPctDJ30indexreturns <- 100*mk.returns(dji,type="relative");
#Create the single Factor which is the dow jones index pct returns
#Through R-2.5.1, timeSeries class originally belong in package fCalendar. 
#Version 221.10065 used cutSeries()method to select data only between the 'to' and 'from' dates. 
#Version 240.10068 used cut(). Both required the day PRIOR to desired start in "from".
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
Fdata <- window(tsRelPctDJ30indexreturns, from="1992-01-01", to="1998-12-30");
Fdata <- as.vector(seriesData(Fdata));


###### Observed Factor Model ####################
#See Table 3.6 and Example 3.36 on pp. 108-9 of QRM book.
#The following function runs a linear regression of the parameter v on the parameter
#Fdata. In regression terms, v is the RESPONSE and Fdata is the predictor. In the factor
#model, v is the large-dimensional variable and FData is reduced-dimensional set of factors.
#Here we are using a single factor (the dji).
# I.e. the individual stock v is the response variable and Fdata (which is the dji
#pct relative return is the PREDICTOR)
vfunc <- function(v,Fdata)
{
#linear regression of the individual stock v [RESPONSE]on the dji [PREDICTOR]Fdata
lm(v ~ Fdata)
}
dim(Xdata)

#The function apply(X, MARGIN, FUN, ...)returns a vector or array or list of values 
#obtained by applying a function FUN passed as 3rd parameter to margins of an array 
#passed as 1st parameter. The 2nd argument MARGIN is 1 if function applied to rows, 
#2 if applied to columns, or c(1,2) if applied to both rows and columns.  The 4th 
#parameter is passed as a named parameter to vfunc().  In other words, the function
#is applied successively to the next column of Xdata each time the vfunc() is called.
# Fit 10 separate univariate regressions: Fdata (dji) vs each individual element of 
#Xdata ("MO", "KO",...)
outa <- apply(Xdata,2,vfunc,Fdata=Fdata);

#Get the names of each result in the list of the ten regressions, e.g. outa$MO will 
#give the name of the 'lm' object returned from lm() for the RESPONSE variable "MO"
names(outa);

# LINEAR REGRESSION MODEL
#The lm() function returns an object of class 'lm'.  The function summary('lm') extracts the
#results for a particular 'lm' object into a class called 'summary.lm'. A 'summary.lm'
#object has the following SLOTS: 
#'$call' gives the regression formula,
#'$residuals' gives the resudual for each position
#'$coefficients' gives the Estimate, StdErrork, t-value, and p-value
#'$sigma' gives the std error
#'$df' gives degrees of freedom
#'$r.squared' gives the r-squared value
#'$adj.r.squared' gives adjusted r-square
#'$fstatistic'
#'$cov.unscaled' 
#which can be accessed via terminology like
#  summary(outa$MO)@r.squared
#  summary(outa$MO)@sigma
 
summary(outa$MO);

#This function runs the regression, fetches the summary() and returns only the
# r-squared value for the regression. In regression, the 'response' is the dependent
# variable y.  The 'predictor' is the (vector of) independent vai
getr2 <- function(response,predictor)
{
    summary(lm(response ~ predictor))$r.squared
}

#return the vector of R-Squares when regressing each element of Xdata against
#the predictor variable 'Fdata' (the dow jones index) separately.  Hence apply()
#causes 10 separate regressions to be run successively, one on each column of Xdata.  
r2 <- apply(Xdata,2,getr2,predictor=Fdata);


# Fit Multivariate Regression.  See p.107 in QRM.  ('Multiple' regression means we have
# multiple x (predictor) variables. 'Multivariate' regression means we have multiple
# y (response) variables. If we have both, "it is multivariate, multiple regression".)
#This example is NOT 'multiple' because there is only one factor in the Fdata
#(the dow jones return).  However, each element in the RESPONSE variable is actually a row
#of returns for each of the ten individual stocks at a point in time. Hence there are 
#effectively 10 response variables so this is the multivariate case. The return class
#is 'mlm' rather than 'lm'. The 'mlm' class is for MULTIPLE RESPONSE MODELS, i.e the
#RESPONSE (y) variable is not just a vector; it is a MATRIX!!!!!  Note we run this ONCE,
#NOT COLUMNWISE as we did previously using apply().
out <- lm(Xdata ~ Fdata);
names(out)
summary(out)
par.ests <- coef(out)
par.ests

# Get parameter estimates
a <- par.ests[1,] #row 1 should be the coefficient associated with the intercept for each RESPONSE.
#If there were multiple factors, B would be a matrix.  Here it is a vector.
B <- par.ests[2,] #row 2 should be the beta coefficients for each RESPONSE.


# Plot betas vs r-squares
plot(B,r2)
text(B,r2,names(r2))
abline(v=1)

#local function
print.loadings.local <- function(x, cutoff=0.1)
{
  cx <- format(round(unclass(x),3))
  nc <- nchar(cx[1])
  cx[abs(x) < cutoff] <- paste(rep(" ",nc), collapse = "")
  print(cx);
}

# Get residuals and estimate covariance matrix for X (XData)from regression:
epsilon <- resid(out);
# Are off-diagonal correlations small ?
round(cor(epsilon),2);
print.loadings.local(cor(epsilon));
# Are the errors uncorrelated with factors (there is only one factor, the dji)
round(cor(epsilon,Fdata),2);

# Construct implied covariance matrix for X (XData). #See equation 3.61 on p. 104 of QRM.
Psi <- var(epsilon) #How close to a diagonal matrix is this? Are errors uncorrelated? See test above.
Omega <- as.matrix(var(Fdata))
#if parameter passed to diag() is MATRIX, then diag() returns a vector containing diagonal
#elements of matrix.  If parameter passed to diag() is VECTOR, then diag() returns a matrix
#with the diagonal elements on the matrix.  In the following, diag(Psi) returns a vector since
#Psi is a matrix.  Hence diag(diag(Psi)) creates a matrix withe the diag(Psi) on the diagonal.
#See equation 3.61 on p. 104 of QRM shows this is covariance matrix for X implied by factor model
#where Psi must be diagonal, i.e. the errors must be uncorrelated with one another 
Sigma <- B %*% Omega %*% t(B) + diag(diag(Psi))
dimnames(Sigma) <- list(dimnames(par.ests)[[2]],dimnames(par.ests)[[2]])
Sigma

# Look at implied correlation matrix.  At the bottom of p. 107 in QRM book, we see written:
# "it is sometimes of interst to form the covariance matrix implied by the factor model and
#compare this with the original sample covariance matrix S of the data."  Here we compare the
#implied correlation matrix to the sample correlation matrix of the data instead.
#R's stat package already has a built-in converter for cov to correlation called cov2cor()
cor.fact <- cov2cor(Sigma);  #CovToCor(Sigma) is S-Plus
cor.data <- cor(Xdata)
round(cor.fact,2)
round(cor.data,2)
print.loadings.local(cor.data - cor.fact)

###################################################################################33
# Principal Components

# Spectral Decomposition
#S-Plus has an unbiased switch. R does not. To set the biased mode, multiply by (n-1/n)
#S <- var(Xdata,unbiased=FALSE)  #S-Plus version
length <- dim(Xdata)[1];
S <- var(Xdata)*(length-1)/length;
rm(length);
#Get the eigenvalues:
tmp.eigen <- eigen(S)
names(tmp.eigen)
L <- diag(tmp.eigen$values)
G <- tmp.eigen$vectors

# Check spectral decomposition
check <- G %*% L %*% t(G)
check[1:3,1:3]
S[1:3,1:3]

# R calculates principal components
out <- princomp(Xdata)
summary(out)
out$loadings; #see the loadings. princomp() returns a list including loadings.

# Check that sum of variances of PCs equals sum of variances of original data
PC.variances <- out$sdev^2
sum.PCvariances <- sum(PC.variances)
sum.PCvariances
traceS <- sum(diag(S))
traceS
cumsum(PC.variances)/sum(PC.variances)

# Graphical displays of explained variances and loadings
screeplot(out);
#The following gives the same output when 'out' is the return value from princomp()
plot(out);
loadings(out)
#Note that the SIGNS of the loadings(out) values are the OPPOSITE of what you get in S-Plus. Why?
print(loadings(out), sort=TRUE)

plot(loadings(out), sort=TRUE)


G <- loadings(out)
FdataTrans <- Xdata %*% G
# could subtract mean returns
# Build two new timeSeries objects using just the first two columns of 
Fact1 <- timeSeries(FdataTrans[,1],seriesPositions(Xdata.ts));
Fact2 <- timeSeries(FdataTrans[,2],seriesPositions(Xdata.ts));
plot(Fact1)
plot(Fact2)


#R does not use Finmetrics modules from S-Plus.  Try to substitute the following
#from the pcaPP library which you can download from CRAN:
library(pcaPP);
outA <-PCAProj(Xdata, k=1);
outA$loadings
outB <-PCAProj(Xdata, k=2); #use two factors
outA$loadings


# Finmetrics contains functions to simplify factor modelling with PCA
#module(finmetrics)
# Model with one factor
#outA <- mfactor(Xdata)
#outA
#loadings(outA)
#plot(outA,which.plots=c(1,2))

# Model with two factors
#outB <- mfactor(Xdata,k=2)
#outB
#loadings(outB)
#plot(outB,which.plots=c(1,2))

# Closer look at factor loadings
#par(mfrow=c(1,2))
#barplot(loadings(outB)[1,],names=colIds(loadings(outB)),horiz=TRUE,main="Factor 1 loadings")
#barplot(loadings(outB)[2,],names=colIds(loadings(outB)),horiz=TRUE,main="Factor 2 loadings")
# Standardize loadings to give portfolio weights
#plot(summary(mimic(outB)))
#par(mfrow=c(1,1))


# Correlation structures implied by factor model
#R.pcafact <- CovToCor(vcov(outA))
#R.pcafact2 <- CovToCor(vcov(outB))
#R <- CovToCor(S)
#round(R.pcafact,2)
#round(R.pcafact2,2)
#round(R,2)

# Correlation structure of residuals/errors
#print.loadings(cor(resid(outA)))
#print.loadings(cor(resid(outB)))
