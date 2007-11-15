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
#Choose only the data from 1/1/1993 to 12/31/2000.  
#Through R-2.5.1, timeSeries class originally belong in package fCalendar. 
#Version 221.10065 used cutSeries()method to select data only between the 'to' and 'from' dates. 
#Version 240.10068 used cut(). Both required the day PRIOR to desired start in "from".
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
partialDJ30daily <- window(partialDJ30dailyTS,from="1993-01-01", to="2000-12-31");

dim(partialDJ30daily@Data);
rm(partialDJ30dailyTS);
partialDJ30dailyMatrix <- seriesData(partialDJ30daily);
partialDJ30weeklyMatrix <- seriesData(aggregateWeeklySeries(partialDJ30daily, FUNC= colSums));
partialDJ30monthlyMatrix <- seriesData(aggregateMonthlySeries(partialDJ30daily, FUNC= colSums));
partialDJ30quarterlyMatrix <- seriesData(aggregateQuarterlySeries(partialDJ30daily, FUNC= colSums));

#########END INITIAL CODE##################

mod.GAUSS <- fit.norm(partialDJ30dailyMatrix )
mod.NIG <- fit.mNH(partialDJ30dailyMatrix ,symmetric=FALSE,case="NIG")
mod.HYP <- fit.mNH(partialDJ30dailyMatrix ,symmetric=FALSE,case="hyp")
mod.t <- fit.mst(partialDJ30dailyMatrix)
mod.NIGs <- fit.mNH(partialDJ30dailyMatrix ,symmetric=TRUE,case="NIG")
mod.HYPs <- fit.mNH(partialDJ30dailyMatrix ,symmetric=TRUE,case="hyp")
round(c(mod.GAUSS$ll.max,mod.t$ll.max,mod.NIGs$ll.max,mod.HYPs$ll.max,mod.NIG$ll.max,mod.HYP$ll.max),1)

######## weekly data
mod.GAUSS <- fit.norm(partialDJ30weeklyMatrix )
mod.NIG <- fit.mNH(partialDJ30weeklyMatrix ,symmetric=FALSE,case="NIG")
mod.HYP <- fit.mNH(partialDJ30weeklyMatrix ,symmetric=FALSE,case="hyp")
mod.t <- fit.mst(partialDJ30weeklyMatrix )
mod.NIGs <- fit.mNH(partialDJ30weeklyMatrix ,symmetric=TRUE,case="NIG")
mod.HYPs <- fit.mNH(partialDJ30weeklyMatrix ,symmetric=TRUE,case="hyp")
round(c(mod.GAUSS$ll.max,mod.t$ll.max,mod.NIGs$ll.max,mod.HYPs$ll.max,mod.NIG$ll.max,mod.HYP$ll.max),1)

