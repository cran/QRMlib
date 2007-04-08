#This is zzz.R.
#This code should be executed immediately if a library(QRMlib) call is executed.
.First.lib <- function(lib, pkg)
{
  library.dynam("QRMlib", pkg, lib,file.ext = .Platform$dynlib.ext);
}

.Last.lib <- function(lib, pkg)
{
  library.dynam.unload("QRMlib", pkg, lib,file.ext = .Platform$dynlib.ext);

}