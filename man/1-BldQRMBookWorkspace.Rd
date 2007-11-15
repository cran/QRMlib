\name{QRMBook-workspace}
\alias{QRMBook-workspace}
\title{How to Build a QRMBook Workspace in R to Use QRMlib}
\description{
Follow these instructions to build a QRMBook workspace in R where you can
run the book's scripts which build most of the graphics plots and tables in
the book\emph{Quantitative Risk Management: Concepts, Techniques and Tools}  
by Alexander  J. McNeil, Rudiger Frey and Paul Embrechts.\cr
\cr

The QRMlib contains scripts which explain topics in most QRM Book chapters.\cr
The folders containing the scripts are named \cr
\dQuote{C:\\Program Files\\R\\R-2.6.0\\library\\QRMlib\\Chap2}, \dQuote{...\\Chap3}, etc. \cr
You may open these scripts from within R by choosing \emph{File | Open Script} from
the R-menu and then moving to the appropriate Chapter script for the QRM Book. Many
chapters contain multiple scripts.
}
\details{
\emph{Instructions to Build the QRMBook workspace}\cr
The following example assumes you are using R version R-2.6.0 in Windows.  
If you are using a different version, substitute your version number in the 
following instructions.\cr

0. Be sure you have R closed.\cr

1. Using MyComputer or Explorer test for the existence of the folder\cr
\emph{C:\\Program Files\\R\\R-2.6.0\\Users}.\cr
If \sQuote{Users} folder does NOT EXIST, create it.\cr 
Each separate project should be built in a subfolder of the \sQuote{Users} folder.\cr
Next create a \sQuote{QRMBook} subfolder beneath the \sQuote{Users} subfolder.\cr 
You should now have a folder\cr
\emph{C:\\Program Files\\R\\R-2.6.0\\Users\\QRMBook}\cr
which you will use only for running code from the QRMlib package.\cr  

2. Right-click the desktop and choose \emph{New | Shortcut} from the menu.\cr

3. Copy the following line (including quotation marks) into your clipboard:\cr 
 \dQuote{C:\\Program Files\\R\\R-2.6.0\\bin\\Rgui.exe}\cr 
and paste the line into the box labeled \dQuote{Type the location of the item}\cr

4. Click the Next> button.\cr

5. Type \emph{QRMBook} (without any quotation marks) into the box labeled \dQuote{Type
a name of this shortcut}.  Then click the Finish button.\cr

6. Find the shortcut you just created on your desktop. It will be 
labeled \dQuote{QRMBook}. Right-click the icon for the shortcut
and choose \sQuote{Properties}.\cr

7. The \sQuote{Start in} box says \dQuote{C:\\Program Files\\R\\R-2.6.0\\bin}.  Modify it to
read \dQuote{C:\\Program Files\\R\\R-2.6.0\\Users\\QRMBook} (be sure to include the 
quotation marks).  Then click OK.
}

\note{
You may now launch the QRMBook workspace by double-clicking the newly-created desktop icon.\cr
This will open R with a workspace pointing to \sQuote{...Users\\QRMBook}.\cr
However,there are still two problems with the workspace:\cr

1. You want to avoid having to load the QRMlib each time you open the workspace. \cr
See \code{\link{profileLoadLibrary}} to resolve this issue\cr

2. You want to use data without issuing the command \emph{data(filename)} each time you open the workspace. \cr
See \code{\link{storeDataInWorkspace}} to resolve this issue.\cr
}

\seealso{
   \code{\link{profileLoadLibrary}} \cr
   \code{\link{storeDataInWorkspace}}
}
\keyword{package}


