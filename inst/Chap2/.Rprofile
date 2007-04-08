#This file uses .first and .Last functions rather than the 
#  options(defaultPackages()) method to set the libraries to be loaded when
# a project starts.  The fCalendar package will also be loaded since it is
# required by QRMlib. Other required libraries include its,Hmisc, chron, mvtnorm
.First <- function()
{
   library(QRMlib)
}

.Last <- function()
{
   detach(package:QRMlib)
}
