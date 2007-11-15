# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# QRMlib 1.4.2
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and danish and DJ data sets##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);

#if you have previously opened the danish (fire loss) and DJ (Dow Jones 30 stocks)
# timeSerie AND aved the workspace, you may comment out the following lines:
data(danish);
data(DJ);
#Alternatively, if you want to load the dataframes instead of timeSeries,
#activate the following lines:
#data(danish.df)
#data(DJ.df);

#################################################


# Analysis of threshold exceedances with GPD
#See example 7.23, p. 280 in QRM
#Losses are expressed in units of 1,000,000 Kroner so they are already over a high threshold.
#However, we will set a still higher threshold
plot(danish, main="Danish fire losses 1980-1990", ylab="Losses in Millions of Kroner",type="l");
grid(); #add gridlines to plot

#See why threshold set to 10 in fit.GPD()? Slight kink below value of 10 (million kroner)?
MEplot(danish, main="Mean Excess Plot (MEP) Danish Fire Losses 1980-1990");
#MEplot(danish[danish>0], main="Mean Excess Plot (MEP) Danish Fire Losses 1980-1990");

losses <- seriesData(danish);

#Run xiplot() to identify how shape parameter changes with other threshold values:
#REQUIRED ARGUMENT: 
# 'data': vector or time series of data 
#OPTIONAL ARGUMENT:
# 'models': number of consecutive gpd models to be fitted 
# 'start': lowest number of exceedances to be considered 
# 'end': maximum number of exceedances to be considered 
# 'reverse': should plot be by increasing threshold (T) or number of extremes (F)
# 'ci': probability for asymptotic confidence band; for no confidence band set to F
# 'auto.scale': whether or not plot should be automatically scaled; if not, xlim and ylim graphical parameters may be entered 
# 'labels': whether or not axes should be labelled 
# 'table': should a table of results be printed?
#For every model "fit.GPD" is called. Evaluation may be slow:
xiplot(danish, models=5, start= 25, end = 200, reverse=TRUE, 
  ci=0.95,auto.scale=TRUE,labels=TRUE,table=TRUE);

#SU: Run hillPlot to show what happens with the Hill Plot.  See Example 7.27, p. 287 in QRM
hillPlot(danish, option = "alpha", start = 5, end = 250, p = 0.99);
hillPlot(danish, option = "alpha", start = 5, end = 60, p = 0.99);


#Fit the GPD in two ways: using MLE and 'pwm"
mod <- fit.GPD(danish,10); #MLE by default
mod$par.ests
#Fit using probability-weighted moments (pwm method)
modb <- fit.GPD(danish,10,method="pwm")
modb$par.ests

#Calculate and display the VaR (quantile) and Expected Shortfall at .99 and .999 for MLE estimators
#We will add these to the tail plots below.  Note the first parameter input is the result of the 
#GPD fit obtained via fit.GPD(danish, 10) so it shows values for losses over 10 kroner:
RMs <- RiskMeasures(mod,c(0.99,0.999))
RMs


#Show the tail plots (a graphic fit comparison):
plotTail(modb, main="Tail Plot for Danish Fire Losses 1980-1990 with 'pwm'-fit and xi=.52");
plotTail(mod, main="Tail Plot for Danish Fire Losses 1980-1990 with 'MLE'-fit and xi=.49");
#Add vertical lines to the MLE tail plot showing the position of the VaR (quantile)
#and the Expected shortfall (sfall):
abline(v=RMs[,"quantile"]) #equivalent to abline(v=RMs[,2]);
#use dotted line (lty=2) for shortfall:
abline(v=RMs[,"sfall"],lty=2) #equivalent to abline(v=RMs[,3],lty=2);


showRM(mod,0.99,RM="VaR");
showRM(mod,0.99,RM="ES");
showRM(mod,0.995,RM="VaR");
showRM(mod,0.995,RM="ES");

#Shape plot for GPD
xiplot(danish);


# Analysis of Microsoft return data
DJreturns <- mk.returns(DJ);
MSFT <- DJreturns[,"MSFT"];
plot(MSFT,type="l");
grid(); #add gridlines to plot

nreturns <- -seriesData(MSFT);
MEplot(nreturns[nreturns>0]);
mod2 <- fit.GPD(nreturns,0.03);
plotTail(mod2);
showRM(mod2,0.995,RM="VaR");
showRM(mod2,0.995,RM="ES");


# Simulated t example
#Use data.dump() in S-Plus to get the S-Plus simulated data here
#data.restore("tsimuldata.dump", print=TRUE)
#tsimuldata <- data.matrix(tsimuldata.df)
#Alternatively, simulate data via rt() here in R. 
tsimuldata <- rt(1000, df = 3)
#Fit to GPD using 100 exceedances rather than a specified threshold:
out <- fit.GPD(tsimuldata, ne = 100)
plotTail(out, extend = 1.5)
showRM(out,0.99,RM="VaR")
abline(v = qt(0.99,3))
showRM(out,0.995,RM="VaR")
abline(v = qt(0.995,3))


#OPTIONAL. Added by SU on 11/20/2006.
#Run these to examine examples 7.24, 7.25, and 7.26 on pp. 280-6 of QRM
#Several examples (7.24, 7.25, and 7.26) in QRM Ch 7 use ATT weekly data.  
#Use ALL complete-week ATT data from 1991-2000.  The first week of the series has only 
#3 days, so we'll discard those and start at 1-6-1991.  This gives us exactly 521
#weeks of data as reported in QRM, p. 281, Example 7.24 (AT&T weekly loss data)
Ret.DJ <- mk.returns(DJ);
#Through R-2.5.1, timeSeries class originally belong in package fCalendar. 
#Version 221.10065 used cutSeries()method to select data only between the 'to' and 'from' dates. 
#Version 240.10068 used cut(). Both required the day PRIOR to desired start in "from".
#R-2.6.0. RMetrics 260.72 moved timeSeries to fSeries from fCalendar. Used window() in place of cut().
#No longer need prior date:
DJ30dailyTSFull <- cut(Ret.DJ, from="1991-01-07", to="2000-12-31");

#Call our function to aggregate daily returns to weekly returns for all DJ stocks:
DJ30weeklyTSFull <- aggregateWeeklySeries(DJ30dailyTSFull);
#Extract the AT&T timeSeries log returns from 1993-2000 log returns for all DJ stocks:
attWeeklyTSFull <- DJ30weeklyTSFull[,"T"];
#Remove the DJ30weeklyTS to save disk space if desired:
rm(DJ30weeklyTSFull);
rm(Ret.DJ);
#Calculate LOSS PERCENTAGES rather than returns:
attWeeklyLossesTSFull <- attWeeklyTSFull;
attWeeklyLossesTSFull@Data <- 100.0*(1.0 - exp(seriesData(attWeeklyTSFull)));
#Use all ATT weekly data from 1/6/1991-12/31/2000 so we get full weeks
length(attWeeklyLossesTSFull@positions); #give length of extracted timeSeries
#Plot the percentage losses.  See figure 7.5a on p. 282.
plot(attWeeklyLossesTSFull,main="ATT Weekly Loss Percentages", ylab="Losses in %",type="l");
grid(); #add gridlines to plot

#Build a Mean Excess Plot.  See figure 7.5b on p. 282.
MEplot(attWeeklyLossesTSFull@Data[attWeeklyLossesTSFull@Data > 0]);
#Choose the threshold at 2.75 from the kink in the previous MEPlot(). Create Figure 7.5c on p. 282
plotFittedGPDvsEmpiricalExcesses(attWeeklyLossesTSFull, threshold=2.75); 
#Create Figure 7.6 on p.285.  This figure contains 3 separate graphs: a plotTail() figure and a
#VaR(0.99) and ES(0.99) graph. The plotTail() and VaR() can be obtained together from a call
#to ShowRM(...,RM="VaR").  The plotTail() and ES() can be obtained together from a call to
#ShowRM(...,RM="ES").  To get all three together, we must call both ShowRM(...)while using 
#split.screen(c(1,11) before the first ShowRM(...,RM="VaR") call.  We then use screen(1,new=FALSE) to indicate
#we do NOT make a new graph but will superimpose the 2nd graph on the 1st.  We then call the 2nd
#ShowRM(...,RM="ES").  Finally we call close.screen(all = TRUE) to indicate we are finished putting
#graphs on top of one another.
mod3 <- fit.GPD(attWeeklyLossesTSFull, threshold=2.75);
split.screen(c(1,1));
showRM(mod3,0.99,RM="VaR");
screen(1,new=FALSE);
showRM(mod3,0.99,RM="ES");
close.screen(all = TRUE); 
#Use the shape plot with a 95% confidence interval to look at shape estimators with
#between 20 and 150 exceedances.  Note that our u=2.75 threshold has 102 exceedances.
#This gives you figure 7.7a on p. 286 of QRM. Be sure to use reverse=FALSE to show graph with
#x as increasing exceedances so we can plot a vertical line at 102:
xiplot(attWeeklyLossesTSFull, start=20, end=150, ci=0.95, auto.scale=TRUE, labels=TRUE, reverse=FALSE);
#Draw vertical line at Nu = 102 exceedances:
abline(v=102, lty=4);
#Draw horizontal line at xi=0.22, the GPD shape parameter fitted from fit.GPD():
abline(h=0.22, lty=2); #or use abline(h=mod3$par.ests[1], lty=2);






