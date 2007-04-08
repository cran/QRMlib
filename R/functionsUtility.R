# QRMlib: version 1.4
# this file is a component of QRMlib 

# Copyright (C) 2006 Alexander McNeil 
# R-language additions Copyright (C) 2006-2007 Scott Ulman

# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software 
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA

# Contact: Alexander J. McNeil:  mcneil@math.ethz.ch 
# R-language contact: Scott Ulman : scottulman@hotmail.com 
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
##########################################################

.on.attach <- function()
{
  cat("QRMlib v1.4, Copyright (C) 2005-2006 Alexander McNeil \n");
  cat("R-language modifications Copyright (C) 2006-2007 Scott Ulman \n");
}

############################################################
#SU 06/15/2006: changed code to work with R instead of S-Plus
# Here is old S-Plus code which uses a data input type of tsdata (for time series
#data.  This data type does NOT EXIST in R.  We must use an alternative data type
# like timeSeries from fCalendar or ts from package stats. 
#mk.returns <- function(tsdata,type="log")
#{
#	times <- positions(tsdata)
#	numbers <- seriesData(tsdata)
#	logr <- function(v)
#	{
#		diff(log(v))
#	}
#      relr <- function(v)
#      {
#            diff(v)/v[1:(length(v)-1)]
#      }
#	returns <- switch(type,
#                          log=apply(numbers, 2, logr),
#                          relative=apply(numbers,2,relr)
#                          )
#	times <- times[-1]
#	timeSeries(returns, times)
#}
# Here is the appropriate R-language code.  You must have loaded the fCalendar
# library via library(fCalendar) prior to calling this method
# tsdata should be a timeSeries class created from the R-Metrics fCalendar, fSeries types
# the default type is 'log'; however you can input either 'log' or 'relative'
mk.returns <- function(tsdata,type="log")
{
      if(is.timeSeries(tsdata)== FALSE) stop("1st parameter not RMetrics-type 'timeSeries' class-also check fCalendar loaded")
      if(type == "log")
         return(returnSeries(tsdata, type="continuous",digits=6))
      else if(type == "relative")
         return(returnSeries(tsdata, type="discrete",digits=6))
      else
         stop("type parameter must be either 'log' or 'relative'")
}
#################################################################

hessb <- function(f, x, ep = 0.0001, ...)
{
  eps <- ep * x
  n <- length(x)
  m <- matrix(0., ncol = n, nrow = n)
  for(i in 1.:n) {
    for(j in 1.:n) {
      x1 <- x
      x1[i] <- x1[i] + eps[i]
      x1[j] <- x1[j] + eps[j]
      x2 <- x
      x2[i] <- x2[i] + eps[i]
      x2[j] <- x2[j] - eps[j]
      x3 <- x
      x3[i] <- x3[i] - eps[i]
      x3[j] <- x3[j] + eps[j]
      x4 <- x
      x4[i] <- x4[i] - eps[i]
      x4[j] <- x4[j] - eps[j]
      m[i, j] <- (f(x1, ...) - f(x2, ...) - f(x3, ...) + f(
                                                           x4, ...))/(4. * eps[i] * eps[j])
    }
  }
  m
}

##########################################################################3
#Empirical (cumulative) distribution function
edf <- function(v, adjust = FALSE)
{
  #Make a copy of the original vector:
  original <- v
  #redo the input vector so it's sorted
  v <- sort(v)
  # Function duplicated() determines which elements of a vector are duplicates of elements
  # with smaller subscripts, and returns a logical vector (a vector with T or F in every position)
  #indicating which elements (rows) are duplicates of other elements.
  # Function cumsum(v) applied to a vector places the cumulative sum of the vector for each respective
  # index. Hence vv= cumsum(v) returns a vector where vv[1] = v[1], vv[2] = v[1]+v[2], vv[3]=v[1]+v[2]+v[3]
  # and so forth. Thus if there have been no duplicates up to index j, the value of vv[j] will be j.
  # If there have been duplicates, the value of vv[j] will be j less the number of duplicates to this point
  # in the sorted array.
  vv <- cumsum(!duplicated(v))
   
  repeats <- tapply(v, v, length)
  add <- rep(cumsum(repeats - 1.), repeats)
  #as.numeric(adjust) either adds 1 if the parameter 'adjust' was T or adds 0 if 'adjust' was F.
  df <- (vv + add)/(length(vv) + as.numeric(adjust))
  #11/20/2006: SU Bug-fix. We must cast the following using as.numeric() or the function returns a 
  #vector with just not the probabilites but also the corresponding quantile and is not readily usable.
  #df[rank(original)]
  as.numeric(df[rank(original)]);
}

#############################################################
#SU: 06/15/2006: This function does not seem to be called from anywhere in the Chapter
# scripts. 'Help' says it may be useful in plotting. 
mk.oldts <- function(timeseries)
{
    #This function converts a 'timeSeries' object to an 'its' object. 'its' stands for
    #irregular time series.
       if(!is.timeSeries(timeseries))
          stop("input must be a timeSeries object");

	data <- seriesData(timeseries)
       #The following three line of S-Plus code must be replaced
	#newtimes <- positions(timeseries)
       #julian <- newtimes@.Data[[1]]
	#oldtimes <- dates(julian, out.format = "d.m.y")
       #by these two lines of R code:
       newtimes <- seriesPositions(timeseries)
       oldtimes <- as.POSIXct(newtimes@Data)

	its(data, oldtimes)
}

###################################################################################

ESnorm <- function(p, mean=0, sd=1)
{
  mean+sd*dnorm(qnorm(p))/(1 - p)
}

#################################################################

qst <- function(p, df, mu=0, sigma=1, scale=FALSE)
{
  quant <- qt(p,df)
  if (scale)
    quant <- quant*sqrt((df - 2)/df)
  mu + sigma*quant
}
  
#############################################

ESst <- function(p, df, mu=0, sigma=1, scale=FALSE)
{
  ES <- (dt(qt(p, df), df)/(1 - p)) * ((df + (qt(p, df))^2)/(df - 1))
  if (scale)
    ES <- ES*sqrt((df - 2)/df)
  mu + sigma*ES
}

#############################################

QQplot <- function(data,position=0.5,reference="normal",...)
  {
    n <- length(data)
    plot.points <- ppoints(n,position)
    func <- switch(reference,normal =qnorm,exp =qexp,student=qt)
    xpoints <- func(plot.points,...)
    plot(xpoints,sort(data),xlab=paste("Theoretical",reference),ylab="Empirical")
    NULL
  }

#############################################################

symmetrize <- function(matrix)
{
	matrix2 <- t(matrix)
	matrix2[lower.tri(matrix2)] <- matrix[lower.tri(matrix)]
	matrix2
}

############################################################################
#The following functions were added by Scott Ulman for R-language only.  

############################################################################

#This function factilitate using an R-language timeSeries object from fSeries to
#characterize 'exceedances', i.e. occurrences above a threshold value).  A prime example 
#would be the highest daily stock prices in a daily timeSeries where evidence indicates 
#clustering occurs.  See p. 117 in QRM.

#This method requires that package 'fCalendar' is loaded to represent the 'timeSeries' object and
#that package 'its' is loaded to represent the irregular times series object which will capture
#the exceedances (e.g. the top 50 observations in a three-year timeSeries).
#require(fCalendar)
#require(its)
#thresholdValue is the value in the timeSeries above which you have 50 'excess' observations. The
#threshold will generally be determined by calling the QRM 'findthreshold()' method.
mk.its.exceedances.tS <- function(timeseries, col=1, thresholdValue)
{ 
  #Extract the data and positions from the timeseries.  Choose the desired column
  data <- seriesData(timeseries);
  nCol <- dim(data)[2];
  if(col < 1 || col > nCol)
     stop("col must be at least 1 and cannot exceed the number of columns of data in timeSeries object");
  #Pick only the column of data input as 2nd parameter:
  data <- data[,col];

  dateIndicators <- seriesPositions(timeseries);
  num <- length(data);
  
  #Create a vector which represents the POSITIONS  within the data series where 
  #the threshold is exceeded.  The maximum length of such a vector would be num
  #but the actual length will be substantially smaller
  exceedance.indices <- (1:num)[data>thresholdValue];
  #Fill exceedance.marks vector with all data items from timeseries which exceed the threshold value:
  exceedance.marks <- data[data>thresholdValue];
  #Use the exceedance.indices to extract dates from seriesPositions(timeseries) when exceedences occurred:
  exceedance.dates <- dateIndicators[exceedance.indices];
  #Convert the exceedance.dates to POSIXct dates:
  exceedance.posix <- as.POSIXct(exceedance.dates);
  #Generate an irregular time series from the values and posixct dates:
  its(exceedance.marks,exceedance.posix);
}

#######################################################################3333
# SU: This version differs from 'mk.its.exceedances.tS' because you pass a data vector without
#any date associations as the 1st parameter.  To associate 'dates' with the vector, you must
#pass a time series as the second parameter.  This series should have the dates you want associated
#with the 'datavector'.  A primary example might be when 'datavector' is generated from a 
#random number generator (RNG) where the parameters were set by fitting a time series.  Then the
#'paralleltimeseriesPos' is the position slot of the timeSeries object.
mk.its.exceedances.vector <- function(datavector, paralleltimeseriesPos, thresholdValue)
{
  #timeseriesPos should be the 'positions' portion of a time series (e.g. rseries@positions),
  #i.e the vector of DATES which the user desires to associate with the input datavector
  # which contains no associated dates.  The datavector may for example, be a set of simulated
  #data built from parameters fitted to an actual time series.  Hence we want to associate the
  #'datavector' with the dates of the associated time series.
  num <- length(paralleltimeseriesPos);

  #Make sure parallelTimeSeries has same length as vector:
  if(length(datavector) != num)
    stop("parallel time series positions must have same length as data input vector");

  #Create a vector which represents the POSITIONS  within the datavector where 
  #the threshold is exceeded.  The maximum length of such a vector would be num
  #but the actual length will be substantially smaller
  exceedance.indices <- (1:num)[datavector>thresholdValue];
  #Fill exceedance.marks vector with all datavector items which exceed the threshold value:
  exceedance.marks <- datavector[datavector>thresholdValue];
  #Use the exceedance.indices to extract dates from paralleltimeseriesPos when exceedences occurred:
  exceedance.dates <- paralleltimeseriesPos[exceedance.indices];
  #Convert the extracted exceedance.dates to POSIXct dates:
  exceedance.posix <- as.POSIXct(exceedance.dates);
  #Generate an irregular time series from the values and posixct dates:
  its(exceedance.marks,exceedance.posix);
}

############################################################
# fCalendar. Version 240.10068 function plot.timeSeries() is not working properly
#when x contains multiple columns of data so we temporarily revert to this version:
#library(its) #required
plot.timeSeriesIts <- function (x, reference.grid = TRUE, lty = 1, ...) 
{
    #The 'its' class separates the data into a matrix and the dates into POSIXt transforms of
    #
    setClass("its", representation("matrix", dates = "POSIXt"))
    #Note this internal function is called only from the end of this file by the code line
    #x.its = .its(x@Data, dates = as.POSIXct(seriesPositions(x)), format = x@format)
    # Hence x is NOT the original timeSeries passed into the function. Instead, it is only the
    #data portion extracted via @Data.  Similarly, dates are POSIXct transformations of the dates
    #extracted from the @positions (or via the seriesPositions(x) method. Hence the dimnames(x)[[1]] 
    #will not really be used. Similarly the format of the date-time variable will be substituted and
    #the default in the .its() method will not actually be used!!!!!
    .its <<- function(x, dates = as.POSIXct(x = strptime(dimnames(x)[[1]], 
        format = "%Y-%m-%d")), names = dimnames(x)[[2]], format = "%Y-%m-%d", 
        ...) {
        if (!is(dates, "POSIXt")) 
            stop("dates should be in POSIX format")
        dates = as.POSIXct(dates)
        if (is.null(dim(x))) {
            dim(x) = c(length(x), 1)
        }
        if (is.null(dimnames(x))) {
            dimnames(x) = list(NULL, NULL)
        }
        if (is.null(dimnames(x)[[1]]) & (nrow(x) > 0)) 
            dimnames(x)[[1]] = 1:nrow(x)
        if (is.null(dimnames(x)[[2]]) & (ncol(x) > 0)) 
            dimnames(x)[[2]] = 1:ncol(x)
        if (!(nrow(x) == length(dates))) {
            stop("dates length must match matrix nrows")
        }
        if (!(ncol(x) == length(names))) {
            stop("names length must match matrix ncols")
        }
        dimnames(x)[[1]] = format(dates, format = format, ...)
        dimnames(x)[[2]] = names
        return(new("its", x, dates = dates))
    }
    .itsPlot <<- function(x, y, colvec = 1:ncol(x), type = "l", 
        ltypvec = 1, lwdvec = 1, yrange, format, at, reference.grid, 
        ...) {
        if (missing(yrange)) {
            ylim = range(x, na.rm = TRUE)
        }
        else {
            ylim = yrange
        }
        firstp = TRUE
        xdates = x@dates
        n = dim(x)[1]
        m = dim(x)[2]
        colveclong = rep(colvec, length.out = m)
        ltypveclong = rep(ltypvec, length.out = m)
        lwdveclong = rep(lwdvec, length.out = m)
        for (i in 1:m) {
            vpoints = c(1, which(!is.na(x[, i])), n)
            xxx = x[, i]
            for (j in 1:ncol(xxx)) {
                if (!firstp) {
                  par(new = TRUE)
                }
                else {
                  firstp = FALSE
                }
                plot(x = xdates[vpoints], y = xxx[vpoints, j], 
                  type = type, col = colveclong[i], ylim = ylim, 
                  lty = ltypveclong[i], lwd = lwdveclong[i], 
                  xaxt = "n", ...)
            }
        }
        if (reference.grid) 
            grid()
        axis.POSIXct(x = xdates[vpoints], side = 1, at = at, 
            format = format)
    }
    `[.its` <<- function(x, i, j, drop, ...) {
        if (match("dates", names(list(...)), 0) > 0) {
            dates = list(...)[["dates"]]
            if (!missing(i)) 
                stop("cannot specify both dates and i")
            if (!is(dates, "POSIXt")) 
                stop("dates should be in POSIX format")
            dates = as.POSIXct(dates)
            i = match(dates, dates(x))
            if (any(is.na(i))) 
                stop("some dates are not found")
        }
        if (missing(drop)) {
            drop = FALSE
        }
        if (missing(i)) {
            i = min(1, nrow(x)):nrow(x)
        }
        if (missing(j)) {
            j = min(1, ncol(x)):ncol(x)
        }
        subx <- x@.Data[i, j, drop = drop]
        dates <- x@dates[i]
        ans <- new("its", subx, dates = dates)
        return(ans)
    }
    x.its = .its(x@Data, dates = as.POSIXct(seriesPositions(x)), 
        format = x@format)
    .itsPlot(x.its, ltypvec = lty, reference.grid = reference.grid, 
        ...)
    invisible(x)
}

#############################################################
#This is the S-Plus version of a function which does NOT exist in R.  Before running
#this function in R, you should have declared the 'signalSeries' class via a call like
#  setClass("signalSeries",
#         representation(data="ANY", positions="numeric", start.position="numeric",
#         end.position="numeric", future.positions="numeric",units="character",
#         title="character", documentation="character", attributes="ANY",
#          units.position="character"));
#If you haven't, this function will do the task for you.

signalSeries <- function(data, positions., units, units.position, from = 1, by = 1)
{
   #PARAMETERS:
      # 'data'      a matrix or other class containing data of ANY type
      #'positions.' is an optional 2nd input parameter. It is NOT the 'positions' slot of the 'signalSeries'
      #             but if input, the slot will be set equal to the values in the input.
      #             Parameter 'positions.', if supplied, overrides from, to, by:
      #IMPORTANT NOTE: if all parameters are missing (ie. the call is y <- signalSeries()), then a new
      #                EMPTY class will be instantiated.

      
      #SU: Added these functions internally.  Although regular functions in S-Plus, they do not exist in R.
      as.rectangular <- function(x)
      {
	if(is.rectangular(x))
		x
	else as.data.frame(x)
      }
      is.rectangular <- function(x)
      {
	(is(x, "character") || is(x, "numeric") || is(x, "complex") || is(x, "factor") ||
		is(x, "matrix") || (is(x, "data.frame") && (!is(x, "data.sheet") ||
		!is.ragged(x))) || (is(x, "array") && (length(dim(x)) == 2)) || (is(
		x, "named") && (is.rectangular(x@.Data))) || is(x, "groupVec") || is(
		x, "seriesVirtual") || is(x, "bdVector") || is(x, "bdFrame"))
      }

      #SU. Since R does not have a 'signalSeries' class, we must insure one has been created in the
      #environment before calling new("signalSeries"). If only a virtual function exists in response
      #to the following getClass() call, we must immediately call setClass():
      CRET <- getClass("signalSeries",.Force=TRUE);
      if(CRET@virtual)
          setClass("signalSeries",
            representation(data="ANY", positions="numeric", start.position="numeric",
            end.position="numeric", future.positions="numeric",units="character",
            title="character", documentation="character", attributes="ANY",
             units.position="character"));
      rm(CRET);
       
	# function to create an empty signalSeries object from 'signalSeries' class 
       #if all parameters are missing.  It returns immediately.
	if(missing(positions.) && missing(data) && missing(units) && 
          missing(units.position) && missing(from) && missing(by)) 
          return(new("signalSeries"))

       #Original S-Plus uses numRows() which doesn't exist in R:
	#if(!missing(positions.) && (length(positions.) != numRows(data)))
       if(!missing(positions.) && (length(positions.) != length(data)))
		stop("Positions and data lengths do not agree")

        #Create a new object of class 'signalSeries': 
	 ret <- new("signalSeries");
        #Original S-Plus code calls asSeriesData(data) which effectively calls asRectangular(data)
        ret@data <- as.rectangular(data);   
 
       #If the 2nd input parameter 'positions.' was not provided, generate positions
	if(missing(positions.)) 
       {
              #S-Plus code needs replacing: no numRows() method in R
		 #len <- numRows(ret@data)
              #ret@positions <- numericSequence(from = from, length = len, by = by)
              len <- length(ret@data);
              ret@positions <- seq(from,from+len*by, by);
	}
	else 
       {
		if(!is(positions., "positionsNumeric"))
			positions. <- as(positions., "numeric")
		ret@positions <- positions.
	}
	if(!missing(units))
		ret@units <- as(units, "character")
	if(!missing(units.position))
		ret@units.position <- as(units.position, "character")
	ret
}

#############################################################
#SU: 10/30/2006: special function to aggregate daily/weekly timeSeries objects in R
# into Quarterly timeSeries objects. Rewritten 3/28/2007 for fCalendar 240.10068.
#**********WARNING***************
#This function has not been extensively tested with daily series which don't start in January
#*********************************
#This reworked function is being included in functionsUtility.R.
#It works with fCalendar version 240.10068 and its somewhat improved functionality
#with timeSequence() functions under R-2.4.1.  Use the old version of
#aggregateQuarterlySeries() if running the old version of fCalendar (221.10065) under R-2.2.1

#This function aggregates the returns from a daily (or weekly) 'timeseries' into quarterly returns.
#  PARAMETERS:
#  	1) timeseries: 
#         a daily or weekly timeSeries from the fCalendar library to be aggregated to quarterly data
#  	2)FUNC 
#         is the function used to aggregate; must be one of the values approved for applySeries(). 
#          Default is 'colSums'.  May also use 'max'
aggregateQuarterlySeries <- function(timeseries,FUNC=colSums)
{
  if(!is.timeSeries(timeseries)) 
       stop("timeseries must be of the timeSeries class from fCalendar package");

  qtrCloseDays <-c(31,30,30,31);  #last days of March, June, Sept, December
  qtrCloseMonths <- c(2,5,8,11); #$mon goes from 0 to 11 so we must use one less than the quarter-ending months as we know them
  qtrStartMonths <- c(0,3,6,9);  #$mon goes from 0 to 11
  calcQtrFromMonth <- function(currentMonth)
  {
    if(currentMonth>= 0 && currentMonth<= 2) 
        qtr <- 1
    else if(currentMonth>= 3 && currentMonth<= 5) 
        qtr <- 2
    else if(currentMonth>= 6 && currentMonth<= 8) 
        qtr <- 3
    else if (currentMonth>= 9 && currentMonth<= 11) 
        qtr <- 4
    else  
       stop("month value in timeSeries not in interval 0-11");
    return(qtr);
  }

  #Figure out the 'to' dates to input into timeSequence(): 
  charvecStartDate <- timeseries@positions[1]; #the first date in the daily timeSeries

  charvecLastDate <- timeseries@positions[length(timeseries@positions)];  #last date in timeSeries
  #print(paste("Start date = ", charvecStartDate));
  #print(paste("End date = ", charvecLastDate));

  #The 'to' series contains the series of ENDING DATES for each quarter.  To build it, we must input
  #the ENDING date of the first quarter as the 'from' parameter and the ENDING date of the last quarter
  #as the 'to parameter in a to <- timeSequence(from=, to=) command.
     # EXAMPLE
     #  > data.frame(from,to)
     #        from         to
     #  1 1994-01-01 1994-03-31
     #  2 1994-04-01 1994-06-30
     #  3 1994-07-01 1994-09-30
     #  4 1994-10-01 1994-12-31
  #Calculate the 'to' and 'from' parameters for the to <-timeSequence()
  fromForTo <- timeLastDayInQuarter(charvecStartDate); 
  toForTo <- timeLastDayInQuarter(charvecLastDate);
  #Build the ending dates for each quarter by using timeSequence(). Note these will need 
  #adjusting since there aren't an equal number of days in each quarter.
  #build to <- timeSequence()
  to <- timeSequence(from=fromForTo, to=toForTo, by="quarter", format="%Y-%m-%d", FinCenter = "GMT");
  #We now have an equally-spaced sequence like this which we need to clean up.  The months are supposed to be those of
  #the quarters {3,6,9,12} and the days are supposed to be the corresponding last days of the month
  #[1] "GMT"
  # [1] [1994-03-31] [1994-07-01] [1994-10-01] [1994-12-31] [1995-03-31] [1995-07-01] [1995-10-01]
  # [8] [1995-12-31] [1996-03-31] [1996-07-01] [1996-10-01] [1996-12-31] [1997-03-31] [1997-07-01]
  #[15] [1997-10-01] [1997-12-31] [1998-03-31] [1998-07-01] [1998-10-01] [1998-12-31] [1999-03-31]
  #[22] [1999-07-01] [1999-10-01] [1999-12-31] [2000-03-31] [2000-07-01] [2000-10-01] [2000-12-31]
  #[29] [2001-03-31] [2001-07-01] [2001-10-01] [2001-12-31] [2002-03-31] [2002-07-01] [2002-10-01]
  #[36] [2002-12-31] [2003-03-31] [2003-07-01] [2003-10-01] [2003-12-31] [2004-03-31]
  quartersLength <- to@Dim; #how many quarters in resulting timeSequence?

  ltest <- strptime(to@Data[1],"%Y-%m-%d");
  startQtr <- calcQtrFromMonth(ltest$mon); #which quarter (1,2,3,4) did we start in?
  

#Start to build the monthVector; it will have 4 items if we started in qtr 1, 3 if in qtr 2, etc.
  monthVector <- switch(startQtr,qtrCloseMonths, c(5,8,11), c(8,11), 11);
  lastDayVector <- switch(startQtr,qtrCloseDays, c(30,30,31), c(30,31), 31);
  remainingQtrs <- quartersLength - length(monthVector);
  fullRepetitions <- as.integer(remainingQtrs/4); #how many times to replicate the qtrCloseMonths vector
  remainder <- remainingQtrs - 4*fullRepetitions; #how many extra quarters do we need to process for a fraction of a year?
  #Replicate the qtrCloseMonths vector fullRepetitions times and combine those to the partial monthVector:
  monthVector <- c(monthVector,rep(qtrCloseMonths,times=fullRepetitions));
  lastDayVector <- c(lastDayVector, rep(qtrCloseDays,times=fullRepetitions));
  #Add the remaining partial elements of qtrCloseMonths:
  if(remainder > 0)
  {
    monthVector <- c(monthVector,qtrCloseMonths[1:remainder])
    lastDayVector <- c(lastDayVector, qtrCloseDays[1:remainder]);
  }
  ltest <- strptime(to@Data[],"%Y-%m-%d");
  #Replace all the months in 'to' with the vector just built so all months will be the last month of a quarter:
  ltest$mon <- monthVector;
  #Replace all the days in 'to' with the vector just built:
  ltest$mday <- lastDayVector;
  #Get rid of time portion in display:
  ltest$hour <- 0;
  #Reset the data:
  to@Data[] <- ltest;


  #Now set the 'from <- timeSequence() values:
  fromForFrom <- timeFirstDayInQuarter(charvecStartDate); 
  toForFrom <- timeFirstDayInQuarter(charvecLastDate);
  #Use the new version of timeSequence() from fCalendar 240.10068 and pass by="quarter":
  from <- timeSequence(from=fromForFrom, to=toForFrom, by="quarter", format="%Y-%m-%d", FinCenter = "GMT");
  #Use the new version of applySeries() from fCalendar 240.10068 and pass by="quarterly":
  applySeries(timeseries,from, to, by="quarterly", FUN=FUNC);
  
}


############################################################################
#SU: 10/30/2006: special function to aggregate daily/weekly timeSeries objects in R
# into Monthly timeSeries objects.  Rewritten 3/28/2007 for fCalendar 240.10068.
#**********WARNING***************
#This function has not been extensively tested with daily series which don't start in January
#*********************************
#This reworked function works with fCalendar version 240.10068 and its somewhat improved functionality
#with timeSequence() functions under R-2.4.1.  Use aggregateMonthlySeries() if running the old version
#of fCalendar (221.10065) under R-2.2.1
#This function aggregates the returns from a daily (or weekly) 'timeseries' into monthly returns.
#  PARAMETERS:
#  	1) timeseries: 
#         a daily or weekly timeSeries from the fCalendar library to be aggregated to monthly data
#  	2)FUNC 
#         is the function used to aggregate; must be one of the values approved for applySeries(). 
#         Default here is 'colSums' rather than the 'colAvg' from the applySeries() method. May also use 'max'.
aggregateMonthlySeries <- function(timeseries,FUNC=colSums)
{
  if(!is.timeSeries(timeseries))
       stop("timeseries must be of the timeSeries class from fCalendar package");
  
  #The new version of timeSequence() in fCalendar 240.10068 now allows us to make series which are on
  #the first day of successive months without all the intervention required in earlier versions. However
  #to make consecutive series on the last days of successive months still takes mighty intervention.
  #We no longer need to correct for the unequal number of days in successive months when using FIRST DAYS.
  #Figure out the 'to' dates to input into timeSequence(): 
  charvecStartDate <- timeseries@positions[1]; #the first date in the daily timeSeries in char form
  charvecLastDate <- timeseries@positions[length(timeseries@positions)];  #last date in timeSeries in char form
  dateFormat="%Y-%m-%d";

  #print(paste("Start date = ", charvecStartDate));
  #print(paste("End date = ", charvecLastDate));
  
  #Generate 'to' and 'from' sequences to pass to timeSequence().  The 'to' sequence should contain the series
  #of ENDING dates for each month.  The 'from' sequence should contain the STARTING DATE for each month.
 
  #The 'to3' series will contain the series of ENDING DATES for each month.  To build it, we must input
  #the ENDING date of the first month as the 'from' parameter and the ENDING date of the last month
  #as the 'to parameter in a to3 <- timeSequence(from=, to=) command.
  #The 'from3' sequence should contain the STARTING DATE for each period.
     # DESIRED EXAMPLE
     #  > 
     #        from3         to3
     #  1 1994-01-01 1994-01-31
     #  2 1994-02-01 1994-02-28
     #  3 1994-03-01 1994-03-31
     #  4 1994-04-01 1994-04-30

  #Using new version of timeSequence() from fCalendar version 240.10068:
  test1 <- as.character(timeFirstDayInMonth(charvecStartDate));
  test2 <- as.character(timeFirstDayInMonth(charvecLastDate));
  #Use the new version of timeSequence() from fCalendar 240.10068 and pass by="month":
  from3 = timeSequence(from = test1, to=test2, by="month", format=dateFormat);


  test3 <- as.character(timeLastDayInMonth(charvecStartDate));
  test4 <- as.character(timeLastDayInMonth(charvecLastDate));
  #Use the new version of timeSequence() from fCalendar 240.10068 and pass by="quarter":
  to3 = timeSequence(from = test3, to=test4, by="month", format=dateFormat);

  #We now have two vectors of equal length, 'to3' and 'from3'. Unfortunately, dates in the'to3' 
  #vector are equally spaced in numbers of days. Since all months do not have the same number
  #of days, the sequence has been created to move ahead into the next month. 
  monthsLength <- to3@Dim; #how many elements in each timeSequence?
  repetitions <- as.integer(monthsLength/12); #how many times to replicate the vector
  remainder <- monthsLength - 12*repetitions; #how many extra months of a year do we need to process for a fraction of a year?
 
  #Vector for 1st day of month:
  firstDayInMonthYr <- c(1,1,1,1,1,1,1,1,1,1,1,1);
  #Note that to@Data$mon is the month in the year indicated by a value of n-1 where n is the way we
  #Americans think about months since months in the fCalendar package run from 0 to 11):
  monthNumbers <- c(0,1,2,3,4,5,6,7,8,9,10,11); 
  lastMonthDayStdYr <- c(31,28,31,30,31,30,31,31,30,31,30,31);

  #The 'to3' series is supposed to hold the LAST DAY of each month; hence we must replace the to3@Data
  #with a sequence containing the last day of each month and the correct month.  
  # The R function timeSequence() uses seq.POXIXt() to build the sequence. The help file says that when
  #using "month", the sequence first advances the month without changing the day: if this results in an 
  #invalid day of the month, it is counted forward into the next month: see the examples. Hence rather than
  #getting February 28 after January 31, we instead get March 03.  We need to fix this.
  #We'll adjust for leap years afterwards.
  #Form a vector of the appropriate length with the last days of the month:
  lastDayVector <- rep(lastMonthDayStdYr, time=repetitions);
  #Now add the remainder (portion of a year to) to the vector:
  if(remainder > 0)
     lastDayVector <-c(lastDayVector,lastMonthDayStdYr[1:remainder]);

  #Some of the MONTHS created from timeSequence are also incorrect for the 'to3' series so we need to
  #replace them.  It's easiest just to replace them all using the 'monthNumbers' vector.
  monthVector <- rep(monthNumbers, times=repetitions);
  #Now add the remainder (portion of a year to) to the vector:
  if(remainder > 0)
    monthVector <-c(monthVector,monthNumbers [1:remainder]);

   
  #Convert the to3@Data to a POSIXlt class which contains the day, month, year, etc.
  ltest <- strptime(to3@Data,"%Y-%m-%d");
  #Reset the day sequence in the vector to coincide to the correct lastDayVector:
  ltest$mday <- lastDayVector;
  #Reset the month sequence to have consecutive months (since the 
  ltest$mon <- monthVector;

  #The 'to3' series is supposed to be the last day of each month. If the year is a leap year
  #and the month is February ($mon == 1 since $mon is n-1 where n is the way normal people specify months)
  #then make the last day of the month the 29th. Must add 1900 to year since a timeDate object has 1900 
  #subtracted.
  #The following can't be vectorized so we need a loop. Got warning messages and failure when trying to run
  #if(leap.year(to3@Data$year+1900) & to3@Data$mon == 1) to3@Data$mday = 29;
  #with vectorizing '&' rather than non-vectorizing '&&'. Both failed.
  for (iter in 1:monthsLength)
  {
    if(leap.year(ltest$year[iter]+1900) && (ltest$mon[iter] == 1))
     ltest$mday[iter] <- 29;
  }

  #Reset the entire array of dates in to3 to end on last day of successive months:
  to3@Data[] <- ltest;
  #data.frame(as.character(from3), as.character(to3));
  
  #Use the new version of applySeries() from fCalendar 240.10068 and pass by="monthly":
  applySeries(timeseries, from3, to3, by="monthly", FUN=FUNC);
}  
############################################################################
#SU: 11/15/2006: special function to aggregate daily timeSeries objects in R
# into weekly timeSeries objects. Rewritten 3/28/2007 for fCalendar 240.10068.
#This function aggregates the returns from a daily 'timeseries' into weekly returns.
#  PARAMETERS:
#  	1) timeseries: 
#         a daily or weekly timeSeries from the fCalendar library to be aggregated to monthly data
#  	2)FUNC 
#         is the function used to aggregate; must be one of the values approved for applySeries(). 
#         Default here is 'colSums' which adds the daily returns to create weekly returns.
aggregateWeeklySeries <- function(timeseries,FUNC=colSums)
{
  if(!is.timeSeries(timeseries))
       stop("timeseries must be of the timeSeries class from fCalendar package");

  #Internal function to insure the input timeSeries is daily so we can aggregate to weekly 
  isDaily <- function(timeseries)
  {
    #Loop through the first five dates.  If the difference between one pair of them is
    #1, then we have a daily time series:
    for (n in 1:4)
    {
     #Get difference in days between successive observations:
     diff <- difftimeDate(timeDate(timeseries@positions[n+1]),
          timeDate(timeseries@positions[n]));
     if(diff == 1) return(TRUE)
     #else next
    }
    return(FALSE)
  }

  #Internal function to find the first Friday in the input series
  pos1stFri <- function(timeseries)
  {
    #browser();
    pos = -1; #return -1 if position not found
    for (name in 1:5)
    {
      #wday is 0=Sun,1=Mon,...6=Sat
      #if wday is FRIDAY, return the corresponding position 
      ltest <- strptime(timeseries@positions[name],"%Y-%m-%d");
      if(ltest$wday == 5) 
      {
        pos = name;
        break;
      }
    }
    #return pos of 1st Friday in daily series
    return(pos);
  }


  #Internal function to determine number of weeks to use:
  nGetInitialWeekCount <- function(timeseries)
  {
    #Determine how many weeks to use for length.out parameter in timeSequence() 
    aveTradingDaysinWeek = 252/52;
 
    #Divide length of daily time series by ave trading days in week and round up if needed
    weeks <- tsLength/aveTradingDaysinWeek;
    nWeeks <- tsLength %/% aveTradingDaysinWeek; 
    if(weeks - nWeeks >= 0.5)
       nWeeks <- nWeeks +1;
    return(nWeeks);    
  }

  if(!isDaily(timeseries))
    stop("timeseries must have daily frequency to do a weekly aggregation with this routine");

  #set the 'from' parameter in timeSequence() as the character data associated with the 
  #first Friday in the daily dataset:
  nFirstFri <- pos1stFri(timeseries);
  if(nFirstFri == -1)
     stop("timeseries has no Friday in 1st five days of timeSeries");
  charvecStart <- timeseries@positions[nFirstFri];
  #Get date format in character form for entire time series:
  charvecFormat <- timeseries@format;
  
  tsLength <- length(timeseries@positions);#Length of daily time series
  
  #Determine how many weeks to use for length.out parameter in timeSequence() 
  nweeks <- nGetInitialWeekCount(timeseries); 

  #create a vector of weekly Friday timeDates (each date a Friday separated by 7 day calendar intervals.
  #input starting Friday, number of weeks and format and specify 'by="week''
  to = timeSequence(from=charvecStart, length.out=nweeks, by="week", format=charvecFormat);
 
  #Get difference in days between last item in 'to' vector and last item
  #in total timeseries successive observations. The 'to' vector contains timeDate objects.
  #This rectifies case where the last item created in timeSequence() is 12/22/2000 and
  #last item in total time series is 12/31/2000 so we could have an additional week ending 12/29/200
  diff <- difftimeDate(timeDate(timeseries@positions[tsLength]),to[to@Dim]);
  if(diff >= 7)
  {
    nweeks <- nweeks + 1;
    #Remake the 'to' timeDate vector with one additional week:
    to = timeSequence(from=charvecStart, length.out=nweeks, by="week", format=charvecFormat);  
  }
  
  #create  a corresponding vector of 'from' values each five days prior to the dates in 'to';
  #these will be dates on Monday. Each element in vectors 'from, to' represents a five-day
  #period from Monday to Friday throughout the years specified.
  #Vectorizing operation: subtract number of seconds in 4 days to get 'from' as a vector series of
  #timeDates preceding dates in 'to' by precisely 5 days. (We are ignoring weekend data since our  
  #'to' series elements are separated by 7 days and each 'from' element is 5 days prior to 'to' element
  from = to - 4*24*3600
  #zippo <- data.frame(as.character(from),as.character(to));
  #Use the new version of applySeries() from fCalendar 240.10068. However, do NOT pass any 'by' parameter.
  #Use 'colSums' to sum the values in each 5-day period as the weekly return.  If you pass from and to, then 
  #'by' is NOT needed. 'by' is needed only for monthly or quarterly series in the new function.
  weeklyTS <- applySeries(timeseries, from, to, FUN=colSums);
  return(weeklyTS);
}

############################################################################

#This is an adaptation of S-Plus aggregateSeries(); it works specifically for signalSeries.
#See aggregateWeeklySeries(), aggregateMonthlySeries(), aggregateQuarterlySeries() for timeSeries.
aggregateSignalSeries <- function(x, pos, AGGFUNC, together = FALSE, drop.empty = TRUE, 
        include.ends = FALSE, adj,offset, colnames, by)
{
     if(class(x) != "signalSeries") 
        stop("x must be of the signalSeries class");
      
     ####### SU: Added these functions internally. #########  
     # Although regular functions in S-Plus, they do not exist in R.
 
      as.rectangular <- function(x)
      {
	if(is.rectangular(x))
		x
	else as.data.frame(x)
      }
      is.rectangular <- function(x)
      {
	(is(x, "character") || is(x, "numeric") || is(x, "complex") || is(x, "factor") ||
		is(x, "matrix") || (is(x, "data.frame") && (!is(x, "data.sheet") ||
		!is.ragged(x))) || (is(x, "array") && (length(dim(x)) == 2)) || (is(
		x, "named") && (is.rectangular(x@.Data))) || is(x, "groupVec") || is(
		x, "seriesVirtual") || is(x, "bdVector") || is(x, "bdFrame"))
      }
      anyMissing <- function(x)
      {
         q <- is.na(x) ;
         for(i in 1:length(q))
         if(q[i]) return(TRUE)
         return(FALSE);
      }
      ############end internal functions added ####################

	# Aggregate series object to new positions
 	#origpos <- positions(x) #S-Plus version
       #R-version must get number of columns via dimension from data slot and then build sequence
       origpos <- seq(1,dim(x@data)[1],1); 
	
	# construct an object of type 'signalSeries' for return
	newobj <- new("signalSeries");
       #Set the new 'signalSeries' slots equal to those from the first parameter passed in
	newobj@units.position <- x@units.position
	newobj@units <- x@units
	newobj@title <- x@title
	newobj@documentation <- x@documentation
	newobj@attributes <- x@attributes
	if(length(origpos) < 1)
		return(newobj)
       
       #build vector holding min and max of origpos:
	rng <- range(origpos);
	# construct positions if 'pos' not passed as parameter
	if(missing(pos)) 
      {
	  pos <- seq(from = rng[1], to = rng[2] + by, by = by);
	  diffpos <- diff(as(pos, "numeric")); #force to 'numeric' type
	}
	else 
      {
	  # if they passed in positions, add another at the end, to make into breakpoints
	  poslen <- length(pos)
         #create lagged differences in positions (defaults to lag by 1)
	  diffpos <- diff(as(pos, "numeric"));
	  if(anyMissing(diffpos) || min(diffpos) < 0)
			stop("Positions for aggregation must be monotonically increasing without NA values")
		if(poslen > 1) 
             {
			dm <- max(diffpos)
			pos[poslen + 1] <- pos[poslen] + dm
			diffpos <- c(diffpos, dm)
		}
		else 
             {
			pos[poslen + 1] <- pos[poslen] + 1
			diffpos <- c(diffpos, 1)
		}
	}
	poslen <- length(pos);
	if(length(pos) < 2)
		return(newobj);

 	# make bins, and take endpoint off positions we are keeping
	posbins <- pos
	pos <- pos[ - length(pos)]

	# get right bin endpoints
	if(include.ends) {
		posbins[1] <- min(posbins[1], rng[1])
		posbins[poslen] <- max(posbins[poslen], rng[2])
	}

       #Bugfix: parameter 'labels=FALSE' must be set to force return of integer codes instead
       #of FACTORS.  Also right=FALSE must be set to get the first interval not to be NA
       #poscut <- cut(origpos, posbins, include.lowest = F, left.includ = T) # S-Plus
	poscut <- cut(origpos, posbins, labels=FALSE, right=FALSE); 

	poscut[is.na(poscut)] <- 0
	if(!drop.empty) {
		poscats <- seq(along = pos)
		whichpos <- T
	}
	else 
	{
             #unique() removes any duplicates; then sort the unique vector and allow only strictly
             #positive values in poscats.  Hence poscats should be sorted, unique copy of poscuts 
		poscats <- sort(unique(poscut))
		poscats <- poscats[poscats > 0]
		whichpos <- poscats
	}
	#Apply the function listed in the 2nd argument to poscats to build new series data.
       #The 2nd parameter defines a function; the 3rd-7th parameters to lapply() are passed
       #to the function specified in the 2nd parameter of lapply() as parameters 2-6.  
	newdat <- lapply(poscats, function(thecat, x, nc, poscut, func, 
			together)
	{
              #CRITICAL:  the function sub()in S-Plus is a subscript operator; in R it is
              #a substitution (grep like) method.  Hence replace it with standard subscript notation.
              cols <- if(together) func(x[ poscut == thecat,  ])
		 #cols <- if(together) func(sub(x, poscut == thecat,  )) 
                else if(nc > 0)
			sapply(1:nc, function(col, x, fun)
			{
				#fun(sub(x,  , col)) #S-Plus version
				fun(x[  , col])
			}
			#, sub(x, poscut == thecat,  ), func) #S-Plus version
			, x[poscut == thecat,  ], func)
			else numeric(0)
		cols <- as.rectangular(cols)
		if(is(cols, "matrix") && numRows(cols) != 1)
			cols <- matrix(cols, nrow = 1)
		cols
	}
	#, x@data, numCols(x), poscut, AGGFUNC, together)
       # dim(x@data)[2] gives the number of columns in a signalSeries object x rather than using 
       #the S-Plus function numCols(x)which does NOT EXIST in R.
       ,x@data, dim(x@data)[2], poscut, AGGFUNC, together) 

	# add offset or adj to positions IF they have been provided as parameters:
	if(!missing(offset)) pos <- pos + offset else if(!missing(adj) && (poslen >
		1)) 
	{
		pos <- pos + adj * diffpos
	}


	# rbind these rows together and put into a series
	newdat <- do.call("rbind", newdat)

	if(missing(colnames))
		#colnames <- colIds(x@data) #S-Plus
		colnames <- x@units
	#if(!is.null(colnames) && (numCols(newdat) == length(colnames))) #S-Plus
	#	colIds(newdat) <- colnames   #S-Plus
       if(!is.null(colnames) && (dim(x@data)[2] == length(colnames)))
		newdat@units <- colnames

	newobj@positions <- pos[whichpos]
	newobj@data <- newdat
	newobj
}
############################################################################

