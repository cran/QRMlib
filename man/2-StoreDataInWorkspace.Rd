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
to load the time series data from the file sp500.rda into the workspace
}
\details{
The scripts in the QRM book use data included with the installation.\cr
The following data files holding \emph{timeSeries} objects are located at\cr
 \emph{C:\\Program Files\\R\\R-2.10.0\\library\\QRMlib\\data}
subfolder. If you examine that folder you may see the data files are compressed in a file
named Rdata.zip.  You may extract the data into the folder if you wish to see the names of
each separate data file by using WinZip or PKZip.  The timeSeries (binary) data files include:
cac40.rda, danish.rda, DJ.rda, dji.rda, ftse100.rda, FXGBP.RAW.rda, hsi.rda, nasdaq.rda, 
nikkei.rda, smi.rda, sp500.rda, spdata.rda, spdata.raw.rda, and xdax.rda.\cr
\cr
In addition to timeSeries objects, you may want to use \emph{dataframe} files.  A set of these
data types has been provided as well.  Any dataframe file can be converted to a timeSeries
using the ConvertDFToTimeSeries() method in functionsUtility.R.\cr
\cr
The dataframe filenames are ASCII readable and are also stored in the data subfolder. They include\cr 
cac40.df.R, danish.df.R, DJ.df.R, dji.df.R, ftse100.df.R, FXGBP.RAW.df.R, hsi.df.R, 
nasdaq.df.R, nikkei.df.R, smi.df.R, sp500.df.R, spdata.df.R, spdata.raw.df.R, and 
xdax.df.R.  Note the dataframe files all have a .R extension meaning they are readable in
ASCII english.  Each dataframe contains a \sQuote{.df.} within its filename.
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


