# S-Plus script developed by Professor Alexander McNeil, mcneil@math.ethz.ch
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.2.1\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#################################################


############# Correlation estimation: a simulation experiment
# Generate 3000 samples of 90-day data from bivariate t distribution
# Calculate correlation estimates in each sample using standard
# sample estimator and estimator based on rank correlation

#  Compare Kendall's tau to standard Pearson correlation coefficient for heavy-tailed
#data simulated from multivariate-t with 3 degrees of freedom and rho = 0.5.
# See fitting the data to true distribution by subdividing the simulated data into
#3000 set of 90 observations each.  For each set, compare the Pearson estimate and the
#Kendall's tau to the known true value of 0.5.  Note the large deviations when using
#standard Pearson's estimator.  

#######SIMPLIFIED APPROACH TO HELP EXPLAIN aggregateSignalSeries() methodology ###########
#To see the logic behind aggregateSignalSeries() when AGGFUNC=pearson, run the following
#simple looping code to see how we generate four different correlation estimators from
#a total simulated bivariate sample with 300 rows and 2 columns (the two different marginals)
sampleTotalLength <- 300;
#Generate bivariate-t sample; it has sampleTotalLength (300) rows and 2 columns:
tdataSim <- rmt(sampleTotalLength,df=3,rho=0.5,d=2);
subsampleLength <- 90;
decimalNumSubsamples <- sampleTotalLength/subsampleLength; 
integerNumSubsamples <- sampleTotalLength%/%subsampleLength;
if(decimalNumSubsamples > integerNumSubsamples) 
   integerNumSubsamples <- integerNumSubsamples + 1; #should be 4
#Create a new vector of length integerNumSubsamples:
pearson.corOffDiagonal <- vector(mode="numeric", length=integerNumSubsamples); 
#Each time through the loop, we increase both the starting index j and ending
#index k by subsample length (e.g. by 90). Hence we are evaluating the cor()
#function for successive samples of size 90 from our total sample.  Hence the 
#first correlation estimate will be on observations 1:90, the 2nd on 91:180,
#the third on 181:270, and the fourth on 271:300 (the last subsample is smaller
#than the first three subsamples).  We will end up with four different correlation
#estimators from our total sample of size 300. These will be output into the
#vector pearson.corOffDiagonal which has length 4.
for(i in 1:integerNumSubsamples)  
{
  j <- subsampleLength *(i-1) + 1;  
  k <- min(subsampleLength *i,sampleTotalLength);
  #run cor() but return only the [1,2] off-diagonal element rather than entire matrix:
  pearson.corOffDiagonal[i] <- cor(tdataSim[j:k,])[1,2]
}
#display the results from dividing the 300 observations into 4 smaller subsets:
pearson.corOffDiagonal; 
################ END SIMPLIFIED APPROACH##############################


#######Use aggregateSignalSeries() to derive a multitude of correlation estimators
#by aggregating a large series into a number of smaller series which can each
#be used to calculate an individual correlation estimator.
set.seed(13)
m <- 90
n <- 3000
#Generate a 'matrix' class of simulated values with 2 columns and m*n rows
dataSim <- rmt(m*n,df=3,rho=0.5,d=2)
#Generate a 'signalSeries' class from the 'matrix'.  This associates a position slot with
#increasing integer values to distinguish the bivariate pairs.
dataSimSS <- signalSeries(dataSim)

#The following function returns only element [1,2] from the correlation matrix for the 
#dataset x (i.e. the off-diagonal correlation in a bivariate model).
pearson <- function(x) cor(x)[1,2]

#This function calculates the kendall tau only for the [1,2] element of a bivariate
#correlation matrix. 
kendall <- function(x){
 sin(pi*cor.test(x[,1],x[,2],method="kendall")$estimate/2)
}

####This MAY take 3 or 4 minutes
pearson.cors <- aggregateSignalSeries(dataSimSS,by=m,together=T,AGGFUNC=pearson);
#Extract the data part only:
pearson.cors.data <- pearson.cors@data;

####This MAY take about 5 minutes
kendall.cors <- aggregateSignalSeries(dataSimSS,by=m,together=T,AGGFUNC=kendall);
kendall.cors.data <- kendall.cors@data;

#Plot the results. See Figure 3.5 and read Example 3.31, p. 98 in QRM book.
par(mfrow=c(2,1))
plot(pearson.cors.data,ylim=c(-1,1),ylab="Pearson")
plot(kendall.cors.data,ylim=c(-1,1),ylab="Kendall-Transform")
par(mfrow=c(1,1))
