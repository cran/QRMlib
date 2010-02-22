#This is zzz.R.
#This code should be executed immediately if a library(QRMlib) call is executed.
.First.lib <- function(libname, pkgname)
{
  library.dynam("QRMlib", pkgname, libname, file.ext = .Platform$dynlib.ext);
}
#The .Last.lib function has been altered from function(lib,pkg) to function(libpath)in R-2.10
.Last.lib <- function(libpath)
{
  #consequently the 2nd (pkg) parameter has been eliminated from library.dynam.unload():
  library.dynam.unload("QRMlib",libpath,file.ext = .Platform$dynlib.ext);
}