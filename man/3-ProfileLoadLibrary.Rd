\name{profileLoadLibrary}
\alias{profileLoadLibrary}
\title{Build .Rprofile File to Load QRM Library in QRMBook Workspace}
\description{
The QRMlib package (QRMlib.dll) must be loaded into your QRMBook workspace before its
functions can be used by scripts.\cr  
The appropriate command to load a package into a workspace is \cr
\emph{library(QRMlib)} \cr
It will be more convenient for you to add a \dQuote{.Rprofile} file to your QRMBook workspace than
to invoke the \emph{library(QRMlib)} command each time you start up.  By adding \emph{.Rprofile} to your
workspace, you will eliminate the need to load the library each time you run the workspace.\cr

See details below for two ways to install a .Rprofile file into your workspace.
}
\details{
The installation program for QRMlib placed an .Rprofile file in the 
\emph{...\\library\\QRMlib\\inst\\Chap2} folder. There is also a useful \emph{README.txt}
in the  \emph{...\\library\\QRMlib} folder.\cr

In the Windows Explorer, merely copy the \emph{.Rprofile} file from the QRM library 
\emph{Chap2} folder into the QRMBook workspace your previously created to run scripts 
(see \code{\link{QRMBook-workspace}}). Once the \emph{.Rprofile} file exists in your 
QRMBook workspace, the QRMlib will automatically be loaded into the workspace.\cr

Alternatively, you can build an \emph{.Rprofile} file in your QRMBook folder using Notepad
or some other text editor.  Just perform the following steps:\cr
0. Close your R-workspace if it is open. \cr

1. Copy the next nine lines of code (starting with .First and ending with \})
into the clipboard.\cr
.First <- function() \cr
\{\cr
   library(QRMlib)\cr
\}\cr
 \cr
.Last <- function()\cr
\{ \cr
   detach(package:QRMlib)\cr
\}\cr
 \cr
2. Open Notepad: left-click the Start button, choose Run and 
type \emph{notepad} into the box. We will try to save a file named 
\dQuote{.Rprofile}.  Note the entire file name is an extension with no 
prefix.  I.e. there will be no letters prior to the \sQuote{.} and the 
type of file is an \dQuote{Rprofile} type spelled with a capital R followed 
by all small letters.\cr
\cr
3. Paste the copied code into Notepad.\cr
\cr
4. In Notepad, choose \emph{File | Save As} from the menu.\cr
\cr
5. In the resulting box, click the \emph{Save as Type} drop-down box 
and choose \emph{All Files}.\cr
(We are NOT saving as a .txt type.)\cr
\cr
6. Paste the path \cr
\dQuote{C:\\Program Files\\R\\R-2.10.1\\users\\QRMBook\\.Rprofile} \cr
into the File name box. Be sure to spell \emph{.Rprofile} 
exactly as shown since R uses case sensitivity in opening files even 
though Windows does not. Use R-2.10.0 or whatever higher version of
R you are using in place of R-2.10.1 if necessary.\cr
\cr
7. Click the \emph{Save} button.\cr
\cr
You may now open your QRMBook workspace and the console should show that the 
QRMlib, fSeries, mvtnorm, chron libraries have been loaded.
}

\note{
When you exit the R program, you will be asked whether to save the current workspace 
environment.  If you choose to do so, the data files which you have opened via
data(filename) calls will be saved in the workspace so you don't need to execute any
subsequent data(filename) calls to get previously-loaded data.  If you do not save the
workspace, you must execute \emph{data(filename)} each time you open the QRMBook workspace.
}
\seealso{
   \code{\link{QRMBook-workspace}} \cr
   \code{\link{storeDataInWorkspace}}
}
\keyword{package}


