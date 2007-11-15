#S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
#R-version adapted by Scott Ulman, scottulman@hotmail.com
# QRMlib 1.4.2
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and spdata.raw data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the spdata.raw data set  AND saved 
#the workspace, you may comment out the following line:
data(spdata.raw);
#if you have previously opened the spdata data set  AND saved 
#the workspace, you may comment out the following line:
data(spdata);

#################################################

#data(spdata.raw);
#Show the data
spdata.raw;
#Attach spdata.raw to search list. This allows you to call the variables
#(like Bobligors) without qualifying them by the database name:
attach(spdata.raw);


# Select S&P B rating data for analysis.
# In the spdata.raw set, Bobligors is the number of issuers with a B credit rating.
#Bdefaults is a vector containing the number of defaults in the group with B credit rating.

# Estimating default probabilities and correlations by simple means.
#momest() returns a VECTOR whose first element is the probability of one firm defaulting;
#the second element is the joint probability of two firms defaulting; the 3rd is the joint
#probability of three firms defaulting;...the tenth element is the joint probability of 
#ten firms defaulting.
momest(Bdefaults,Bobligors);

#get probability of single default:
pi.B <- momest(Bdefaults, Bobligors)[1];
#get joint probability of two defaults:
pi2.B <- momest(Bdefaults, Bobligors)[2];
rhoY.B <- (pi2.B-pi.B^2)/(pi.B-pi.B^2);

# Fitting parametric distributions by MLE
mod0 <- fit.binomial(Bdefaults, Bobligors);
mod0;
mod1 <- fit.binomialBeta(Bdefaults, Bobligors);
mod1;

# A little patience is require for the next two models ...

mod2 <- fit.binomialProbitnorm(Bdefaults, Bobligors);
mod2;
mod3 <- fit.binomialLogitnorm(Bdefaults, Bobligors);
mod3;
c(mod0$maxloglik, mod1$maxloglik, mod2$maxloglik, mod3$maxloglik);


# Use Splus correlated data library.  This library IS NOT AVAILABLE IN R.  Substitute an R-methodology.
#library(correlatedData)
#options(contrasts=c("contr.treatment","contr.poly"))
#results <- glme(cbind(defaults,firms-defaults) ~ -1 + rating, random = ~1| year, 
#     family=binomial(probit), data=spdata)

# R-methodology
library(MASS);
library(nlme);
#show the data being used: It is spdata rather than spdata.raw.
spdata;
#the data being used is spdata which has 100 rows and 4 columns: 
#'year', 'rating', 'firms', defaults'
#Use R- library MASS to get glmmPQL
#'year' -'ratings' determine the unique results(20 years 1981-2000 
# with 5 obligor class ratings each year
results <- glmmPQL(cbind(defaults,firms-defaults) ~ -1 + rating, 
  random = ~1| year, family=binomial(probit), data=spdata);
summary(results);
summary(results)$tTable[,1];

detach("package:nlme");
detach("package:MASS");



