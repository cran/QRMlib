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
#partialDJ30daily <- cutSeries(partialDJ30dailyTS,from="1992-12-31", to="2000-12-31");
partialDJ30daily <- cut(partialDJ30dailyTS,from="1992-12-31", to="2000-12-31");

dim(partialDJ30daily@Data);
rm(partialDJ30dailyTS);

#Now create the aggregate weekly, monthly, and quarterly values by summing the daily values
partialDJ30weeklyTS <- aggregateWeeklySeries(partialDJ30daily, FUNC= colSums);
partialDJ30monthlyTS <- aggregateMonthlySeries(partialDJ30daily, FUNC= colSums);
partialDJ30quarterlyTS <- aggregateQuarterlySeries(partialDJ30daily, FUNC= colSums);

# *****EXTRACT MATRIX FROM TIME SERIES, effectively DISCARDING TIME POSITIONS****
partialDJ30dailyMatrix <- seriesData(partialDJ30daily);
partialDJ30weeklyMatrix <- seriesData(partialDJ30weeklyTS);
partialDJ30monthlyMatrix <- seriesData(partialDJ30monthlyTS);
partialDJ30quarterlyMatrix <- seriesData(partialDJ30quarterlyTS);
dim(partialDJ30dailyMatrix);
dim(partialDJ30weeklyMatrix);
dim(partialDJ30monthlyMatrix);
dim(partialDJ30quarterlyMatrix);

### Calculating sample moments. Apply() the mean() function across columns (2nd parameter=2)
Xbar <- apply(partialDJ30monthlyMatrix,2,mean)
Xbar

#The following S-Plus code calculates both biased and unbiased versions of variance:
#S <- var(partialDJ30monthlyMatrix) #,unbiased=F)
#Su <- var(partialDJ30monthlyMatrix,unbiased=T)
#S
#R-language has no switch for biased version: only unbiased version produced:
Su <- var(partialDJ30monthlyMatrix);


R <- cor(partialDJ30monthlyMatrix);
R;

#Calculate the Mahalonobis distance:
# Function returns the squared Mahalanobis distance of all rows in 'x' and
# the vector mu='center' with respect to Sigma='cov'. This is (for vector 'x') defined as
#                 D^2 = (x - mu)' Sigma^{-1} (x - mu)
D <- mahalanobis(partialDJ30monthlyMatrix,Xbar,Su,inverted=F);
#Plot a histogram of the mahalonobis distance:
hist(D);
