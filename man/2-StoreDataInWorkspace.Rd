\name{storeDataInWorkspace}
\alias{storeDataInWorkspace}
\title{How to Store Data in a QRMBook Workspace}
\description{
Data files must be loaded into your workspace before they can be used by scripts.\cr  
The appropriate command to load data into a workspace is \cr
\emph{data(filename)} \cr
where filename is the name of one of the data files WITHOUT its R extension.\cr  
Hence use \cr
\emph{data(sp500)}\cr
to load the data from the file sp500.R into the workspace
}
\details{
The scripts in the QRM book use data included with the installation.\cr
The following data files are located at\cr
 \emph{C:\\Program Files\\R\\R-2.4.1\\library\\QRMlib\\data}
subfolder. If you examine that folder you may see the data files are compressed in a file
named Rdata.zip.  You may extract the data into the folder if you wish to see the names of
each separate data file by using WinZip or PKZip.  The data files include:
cac40.R, danish.R, DJ.R, dji.R, ftse100.R, FXGBP.RAW.R, hsi.R, nasdaq.R, nikkei.R, smi.R,
sp500.R, spdata.R, spdata.raw.R, and xdax.R.
}

\note{
When you exit the R program, you will be asked whether to save the current workspace 
environment.  If you choose to do so, the data files which you have opened via
\emph{data(filename)} calls will be saved in the workspace so you don't need to execute any
subsequent \emph{data(filename)} calls to get previously-loaded data.  If you do not save the
workspace, you must execute \emph{data(filename)} each time you open a script in QRMBook workspace.
}

\seealso{
   \code{\link{profileLoadLibrary}}\cr
   \code{\link{QRMBook-workspace}}

}
\keyword{package}


