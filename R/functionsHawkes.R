# QRMlib: version 1.4.2
# this file is a component of QRMlib 

# Copyright (C) 2005-06 Alexander McNeil 
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

# Contact: Alexander J. McNeil:  a.j.mcneil@hw.ac.uk
# R-language contact: Scott Ulman : scottulman@hotmail.com 
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
#########################################################
#SU.  11/7/2006. R-language modifications require substantial rewrite
#R does not contain a "series" type data item.  Hence data parameter must represent
#either vector, a single column from a matrix or data.frame, or a timeSeries with
#a single data column. 
#This method returns an MPP class which contains the following list:
# list(times=times, marks=marks, starttime=starttime, endtime=endtime, threshold=threshold)
extremalPP <- function(data,threshold=NA, nextremes=NA)
{
    #In R, the data parameter is either a timeSeries with a single item in the data slot
    # or a vector of data; there is no "series" object.  If you pass in a matrix or data.frame
    # the program will select column 2 since column 1 is generally the DATE item.
    #object as in S-Plus
    #if(is(data, "series")) 
    if(is.timeSeries(data)) 
    {
      isTS <- TRUE;
      pos <- seriesPositions(data);
      #SU: In the S-Plus code, 'alltimes' represents the julian day count 
      #between the date in the vector and the origin = '1/1/1960'.  R-language julian 
      #day counts are relative to origin='1/1/1970'. Hence we must ADD 3653 days to the 
      #R-day count to replicate S-Plus results. Replace the following S-Plus line:
      #alltimes <- pos@.Data[[1]]
      #by the next three lines to get an integer vector.  
      tD <- timeDate(data@positions)
      #Get the julian counters for the number of days since origin (S-Plus used by default)
      alltimes <- julian.timeDate(tD) + 3653 #this is a difftime object
      #convert to integer vector with cast. The integer cast is necessary to get rid of the
      # comment prepending the integer which says "the number of days in the timeDiff is".
      alltimes <- as.integer(alltimes);
      #11/8/2007: we must convert "values" to a vector so replace the following line with 
      #its successor:
      #values <- seriesData(data);
      values <- as.vector(data@Data[ ,1]); #using first data column if there is more than one
      #startime is first julian day count reduced by 1:
      starttime <- alltimes[1]-1;
      endtime <- alltimes[length(alltimes)];
    }
    #or it is an equally-spaced series of numeric values:
    else
    {
      isTS <- FALSE;
      #11/08/2007: added cases for vector, data.frame, matrix
      if(is.vector(data))
          values <- data
      else if (is.data.frame(data))
          values <- data[ , 2] #using 2nd column since DATE is usually column 1
      else if (is.matrix(data))
          values <- data[ , 1]
      else
           stop("Input data is not appropriate vector, data.frame, or matrix");
      n <- length(values);
      alltimes <- (1:n)/n;
      starttime <- 0;
      endtime <- 1;
    }
    if(is.na(nextremes) && is.na(threshold))
      stop("Either a threshold or the number of upper extremes must be supplied.")
    if(!is.na(nextremes)){
      threshold <- findthreshold(values, nextremes);}
    #Note the 'times' represent the julian day count since 1/1/1960 (S-Plus origin) when
    #each threshold exceedance occurs when inputting a timeSeries; they represent equally spaced
    #times between 0 and 1 when the input data is merely a vector:
    times <- alltimes[values>threshold]
    #Note the 'marks' represent the value of the exceedance, i.e. the amount by which the
    #value exceeds the threshold 
    marks <- values[values>threshold]-threshold
    out <- list(times=times, marks=marks, starttime=starttime, endtime=endtime, threshold=threshold)
    #Indicate that this class will inherit from the MPP (marked-point-process) class
    oldClass(out) <- c("MPP")
    out
}
########################################################
#This function takes a Marked Point Process (MPP) item, removed the 'marks' slot, and changes
#the name of the class to "PP" from "MPP"
unmark <- function(PP){
  if (class(PP) != "MPP") stop("Not a marked point process");
  PP$marks <- NULL;
  oldClass(PP) <- "PP";
  PP;
}

############################################################
#SU: Due to R warning about S3 generic method consistency, we rename initial parameter as x
#consistent with the plot(x) routine. Also, we pass ... for additional arguments to generic plot().
#According to R-help on plot, x may be the coordinates of points in the plot OR a single plotting 
#structure or function.  Note the value PP must be
#a Marked Point Process
plot.MPP <- function(x, ...){
  if (class(x) != "MPP") stop("Not a marked point process");
  starttime <- x$starttime;
  endtime <- x$endtime;
  times <- x$times;
  marks <- x$marks;
  #SU 11/7/2006: no need to use irregular time series:
  #PPits <- its(marks,times)
  #ts.plot(PPits,xlim=range(starttime,endtime),type="h",ylab="MPP")
  #plot(times, marks, xlime=c(starttime,endtime), ylab="MPP");
  plot(times, marks, xlim=c(starttime,endtime), ylab="MPP", type="h", ...);
}

##############################################################
# Plot an UNMARKED POINT PROCESS
#SU: Due to R warning about S3 generic method consistency, we rename initial parameter as x
#consistent with the plot(x) routine. Also, we pass ... for additional arguments to generic plot().
#According to R-help on plot, x may be the coordinates of points in the plot OR a single plotting 
#structure or function.
plot.PP <- function(x,...)
#plot.PP <- function(PP, ...)
{
  if (class(x) != "PP") stop("Not a point process. Did you call unmark()?");
  starttime <- x$starttime
  endtime <- x$endtime
  #times when exceedances occurred.  We have removed marks which gave size of exceedances.
  times <- x$times  
  #In S-Plus, stepfun() is stepfun(datax, datay, type="left"). The parameters 'datax' 
  #and 'datay' may be of the same length and give the locations of the jumps. Values 
  #of 'datax' MUST BE in INCREASING ORDER. They are here since datax=times will be the 
  #Julian data counts on successive times exceedances occurred.
  #However, in R, 'datay' must be a vector one longer than x (so there is a y[0]).  
  #Then we want x[1] as the value between y[1]-y[0], x[2] as the value between y[2]-y[1], etc. 
  #Hence we build a vector of values increasing by 1 like a count:1,2,3...
  augLength <- length(times)+1;
  count <- 1:augLength;
  #SU: 11/7/2006. This function is currently deprecated in S-Plus. It has a different syntax
  # in R from S-Plus where it has this syntax: stepfunc <- stepfun(times,count,type="left")
  #'times' should be numeric vector giving the knots or jump locations of the step function for stepfun()
  #'count' should be a numeric vector one longer than x, giving the heights of the function values between the x values
  #which effectively means the number of exceedances at that time since the count increases by 1 each time an exceedance
  #occurs
  stepfunc <- stepfun(times,count, f = 0); 
  plot.stepfun(stepfunc,xlim=range(starttime,endtime),ylim=range(0,augLength),
      #use pch='' to indicate we want NO character for drawing the individual points; otherwise we  get little circles.
      xlab="Time",ylab="N events", pch='', main="Time Pattern of Exceedances",...);
  NULL
}

############################################################

fit.POT <- function(PP,markdens="GPD"){
  starttime <- PP$starttime
  endtime <- PP$endtime
  times <- PP$times
  marks <- PP$marks
  span <- endtime-starttime
  par.ests <- length(times)/span
  names(par.ests) <- c("lambda")
  par.ses <- sqrt(par.ests/span)
  names(par.ses) <- names(par.ests)
  names
  ll.max <- -par.ests*span + length(times)*log(par.ests)
  if (class(PP)=="MPP")
    {
    mark.model <- switch(markdens,GPD=fit.GPD(marks,0))
    par.ests <- c(par.ests,mark.model$par.ests)
    names(par.ests) <- c("lambda",names(mark.model$par.ests))
    par.ses <- c(par.ses,mark.model$par.ses)
    names(par.ses) <- names(par.ests)
    ll.max <- ll.max+mark.model$ll.max
  }  
  list(par.ests=par.ests,par.ses=par.ses,ll.max=ll.max)    
}

###############################################################################

sePP.negloglik <- function(theta, PP, case) 
 {
   theta <- abs(theta)
   times <- PP$times
   marks <- PP$marks
   if (class(PP) != "MPP") marks <- 0.0
   endtime <- PP$endtime
   starttime <- PP$starttime
   mu <- theta[1]
   phi <- theta[2]
   voltheta <- theta[-c(1,2)]
   evol <- volfunction(times,times,marks,voltheta,case)
   lambda.contrib <- mu+phi*evol
   term1 <- sum(log(lambda.contrib))
   gamma <- theta[3]
   delta <- switch(case,0.0,theta[4],0.0,theta[5])
   rho <- switch(case,0.0,0.0,theta[4],theta[4])
   if (case <=2)
     terminsum <- (1-exp(-gamma*(endtime-times)))/gamma
   else
     terminsum <- gamma*(1-(1+(endtime-times)/gamma)^(-rho))/rho
   terminsum <- (1+delta*marks)*terminsum
   term2 <- mu*(endtime-starttime)+phi*sum(terminsum)
   out <- term2-term1
   out 
 }

#####################################################################################

seMPP.negloglik <- function(theta, PP, case, markdens) 
 {
   theta <- abs(theta)
   times <- PP$times
   marks <- PP$marks
   endtime <- PP$endtime
   starttime <- PP$starttime
   mu <- theta[1]
   phi <- theta[2]
   voltheta <- theta[-c(1,2,(length(theta)-2),(length(theta)-1),length(theta))]
   evol <- volfunction(times,times,marks,voltheta,case)
   lambda.contrib <- mu+phi*evol
   term1 <- sum(log(lambda.contrib))
   xi <- theta[length(theta)-2]
   beta <- theta[length(theta)-1]
   alpha <- theta[length(theta)]
   scale <- beta + alpha*evol
   markdensfunc <- switch(markdens,GPD=dGPD)
   lambda.x <- markdensfunc(marks,xi,scale,logvalue=TRUE)
   term3 <- sum(lambda.x)
   gamma <- theta[3]
   delta <- switch(case,0.0,theta[4],0.0,theta[5])
   rho <- switch(case,0.0,0.0,theta[4],theta[4])
   if (case <=2)
     terminsum <- (1-exp(-gamma*(endtime-times)))/gamma
   else
     terminsum <- gamma*(1-(1+(endtime-times)/gamma)^(-rho))/rho
   terminsum <- (1+delta*marks)*terminsum
   term2 <- mu*(endtime-starttime)+phi*sum(terminsum)
   out <- term2-term1-term3
   out 
 }

####################################################################

volfunction <- function(anytimes,times,marks,theta,model){  
#Changes by Scott Ulman on 06/02/2006 to convert to R rather than S-Plus. 
#Function called via .C() in QRMsep.c is
# void SEprocExciteFunc(long*n, double *times, long *nmarks, double *marktimes, 
# double *marks, double *beta, long *model, double *result)  
   # Get length of anytimes variable 
   nVecLength = length(anytimes)
   #Get length of times variable
   nevents = length(times)
   #fill a vector of length veclength with 0s
   vector = rep(0,nVecLength)
   tempobj <- .C("SEprocExciteFunc",
              as.integer(nVecLength),
              as.double(anytimes),
              as.integer(nevents),
              as.double(times),
              as.double(marks),
              as.double(theta),
              as.integer(model),
              result = as.double(vector),
              PACKAGE="QRMlib") $result
    #tempobj <- .C("SEprocExciteFunc",
    #             n=length(anytimes),
    #             anytimes=anytimes,
    #             nevents=length(times),
    #             times=times,
    #             marks=marks,
    #             theta=theta,
    #             model=model,
    #             result=rep(0,length(anytimes)),
    #             CLASSES=c("integer","double","integer","double","double","double","integer","double"))
    #             tempobj$result
    #The following will return only '$result' as 'tempobj'
    return(tempobj)
  }

#####################################################################################
#SU: Due to R warning about S3 generic method consistency, we rename initial parameter as x
#consistent with the plot(x) routine. According to R-help on plot, x may be the coordinates
#of points in the plot OR a single plotting structure or function.
#Note that a self-exciting point process created by fit.sePP is a list with elements
#out <- list(PP=PP,par.ests=par.ests,par.ses=par.ses,tstat=tstat,model=model,mark.model=FALSE,
#mark.influence=mark.influence,case=case,converged=fit$message,ll.max=ll.max)
#plot.sePP <- function(fit, ...)
plot.sePP <- function(x, ...)
{
  #PP <- fit$PP
  PP <- x$PP;
  starttime <- PP$starttime
  endtime <- PP$endtime
  times <- PP$times
  marks <- PP$marks
  #theta <- fit$par.ests
  #case <- fit$case
  theta <- x$par.ests
  case <- x$case
  anytimes <- starttime:endtime
  voltheta <- theta[-c(1,2)]
  evol <- volfunction(anytimes,times,marks,voltheta,case)
  intensity <- theta[1]+theta[2]*evol
  plot(anytimes,intensity,type="l",xlim=range(starttime,endtime),xlab="Time",ylab="Intensity",...)
  abline(h=length(times)/(endtime-starttime))
}

#############################################################
fit.sePP <- function(PP,model="Hawkes",mark.influence=TRUE,std.errs=FALSE)
{
  starttime <- PP$starttime
  endtime <- PP$endtime
  times <- PP$times
  marks <- PP$marks
  span <- endtime-starttime
  rate <- length(times)/span
  if (!((class(PP) =="MPP") & mark.influence))
    {
    case <- switch(model,Hawkes=1,ETAS=3)
    mark.influence=FALSE
  }
  else
    case <- switch(model,Hawkes=2,ETAS=4)
  
  #This is the first parameter to pass to nlminb(). Depending on 'case', it has 3-5 parameters:
  theta <- switch(case,c(rate,0,0.1), c(rate,0,0.1,0), c(rate,0,0.1,0.1),c(rate,0,0.1,0.1,0))

  #The R-version of nlminb() has a LIST of control parameters 'control = list()' in lieu of nlminb.control().
  #The R defaults for eval.max and iter.max are precisely those given in the following call, so we can omit
  #sending the control list. PP is an optional parameter which can be set to MPP or PP. 'case' tells whether
  #to pass 'Hawkes' or 'ETAS' for solution. Comment out S-Plus version and substitute R-language version of call:
  #fit <- nlminb(theta.start, sePP.negloglik, control=nlminb.control(eval.max=200,iter.max=150), PP=PP,case=case) 
  #Notice the first parameter passed to nlminb() must be the first parameter of the objective function; the other
  #parameters of the objective function are passed as PP and case:  
  fit <- nlminb(start=theta, objective=sePP.negloglik, PP=PP,case=case);

  #The return values for R are called fit$par so replace the following line with an R-version:
  #par.ests <- fit$parameters
  par.ests <- fit$par;
  par.ests <- abs(par.ests)
  ll.max <- -sePP.negloglik(par.ests,PP,case)
  nms <- c("mu","phi")
  addnms <- switch(case, c("gamma"),c("gamma","delta"),c("gamma","rho"),c("gamma","rho","delta"))
  names(par.ests) <- c(nms,addnms)

  if (std.errs){
    #hessb() is in functionsUtility.R if needed.
    hessian <- hessb(sePP.negloglik,par.ests, PP=PP,case=case)
    varcov <- solve(hessian)
    par.ses <- sqrt(diag(varcov))
    names(par.ses) <- names(par.ests)
  }
  else{
    par.ses <- rep(NA,length(par.ests))
    names(par.ses) <- names(par.ests)
  }
  tstat <- par.ests/par.ses
  out <- list(PP=PP,par.ests=par.ests,par.ses=par.ses,tstat=tstat,model=model,mark.model=FALSE,mark.influence=mark.influence,case=case,converged=fit$message,ll.max=ll.max)
  oldClass(out) <- "sePP"
  out
}

##################################################################################

fit.seMPP <- function(PP,markdens="GPD",model="Hawkes",mark.influence=TRUE,predictable=FALSE,std.errs=FALSE)
{
  if (class(PP) != "MPP") stop("Not marked point process data")
  marks <- PP$marks
  groundmod <- fit.sePP(PP,model,mark.influence,std.errs)
  par.ests <- groundmod$par.ests
  nms <- names(par.ests)
  par.ses <- groundmod$par.ses
  tstat <- groundmod$tstat
  ll.max <- groundmod$ll.max
  converged <- groundmod$converged
  case <- groundmod$case
  mark.model <- switch(markdens,GPD=fit.GPD(marks,0))
  par.ests <- c(par.ests,mark.model$par.ests)
  names(par.ests) <- c(nms,names(mark.model$par.ests))
  if (std.errs)
    par.ses <- c(par.ses,mark.model$par.ses)
  ll.max <- ll.max+mark.model$ll.max
  if (!(mark.model$converged))
    converged <- FALSE
  if (predictable)
    {
    nms <- names(par.ests)
    #The R-version of nlminb() has a LIST of control parameters 'control = list()' in lieu of nlminb.control().
    #The R defaults for eval.max and iter.max are precisely those given in the following call, so we can omit
    #sending the control list. PP is an optional parameter which can be set to MPP or PP. 'case' tells whether
    #to pass 'Hawkes' or 'ETAS' for solution. Comment out S-Plus version and substitute R-language version of call:
    #theta.start <- c(par.ests,0)
    #fit <- nlminb(theta.start, seMPP.negloglik, control=nlminb.control(eval.max=200,iter.max=150), PP=PP,case=case,markdens=markdens)
    theta <- c(par.ests,0); 
    fit <- nlminb(start=theta, objective=seMPP.negloglik, PP=PP,case=case, markdens=markdens);
    #The return values for R are called fit$par so replace the following line with an R-version:
    #par.ests <- fit$parameters
    par.ests <- fit$par
    par.ests <- abs(par.ests)
    names(par.ests) <- c(nms,"alpha")
    ll.max <- -seMPP.negloglik(par.ests,PP,case,markdens=markdens)
    converged <- fit$message
    if (std.errs){
       hessian <- hessb(seMPP.negloglik,par.ests, PP=PP,case=case,markdens=markdens)
       varcov <- solve(hessian)
       par.ses <- sqrt(diag(varcov))
    }
  }
  if (!std.errs)
    par.ses <- rep(NA,length(par.ests))
  names(par.ses) <- names(par.ests)
  tstat <- par.ests/par.ses
  out <- list(PP=PP,par.ests=par.ests,par.ses=par.ses,tstat=tstat,model=model,
       mark.model=TRUE,mark.influence=mark.influence,case=case,predictable=predictable,
        converged=converged,ll.max=ll.max);
  oldClass(out) <- "sePP"
  out
}

####################################################################################

stationary.sePP <- function(sePP)
  {
    mark.model <- sePP$mark.model
    mark.influence <- sePP$mark.influence
    if (mark.model)
      if (sePP$predictable)
        stop("Only implemented for unpredictable marked point processes")
    not.mark.model <- !(mark.model)
    if (not.mark.model & mark.influence) stop("Need auxiliary model for marks")
    par.ests <- sePP$par.ests
    phi <- par.ests[2]
    gamma <- par.ests[3]
    case <- sePP$case
    if(sePP$model=="Hawkes")
      eta <- phi/gamma
    else{
      rho <- switch(case,0.0,0.0,par.ests[4],par.ests[4])
      eta <- phi*gamma/rho
    }
    if (mark.influence) {
      delta <- switch(case,0.0,par.ests[4],0.0,par.ests[5])
      xi <- par.ests[length(par.ests)-1]
      beta <- par.ests[length(par.ests)]
      eta <- eta*(1+delta*beta/(1-xi))
    }
    cluster.size <- 1/(1-eta)
    out <- c((eta<1),eta,cluster.size)
    names(out) <- c("stationary","eta","cluster.size")
    out
  }
