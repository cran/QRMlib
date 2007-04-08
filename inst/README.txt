QRMlib library: version 1.4
Copyright (C) 2005-2006 Alexander McNeil
R-language modifications copyright (C) 2006-2007 by Scott Ulman.

This program contains R and C routines to accompany the book 
Quantitative Risk Management: Concepts, Techniques and Tools.
by McNeil, frey and Embrechts (PUP, 2005)
The C code uses routines in the GSL library.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Contact: Scott Ulman:  scottulman@hotmail.com

## Windows Instructions for R users:
The following instructions pertain to a user with version 2.4.1 or higher 
of R installed.  
************
**** Creating a QRMBook Workspace (User) ***********
1. Using MyComputer or Windows Explorer, create a QRMBook subfolder beneath 
	C:\Program Files\R\R-2.4.1\Users. 
(wherever your R Users folder resides).
2. Right-click the desktop and choose New | Shortcut from the menu.
3. Copy the following line (including quotation marks) into your clipboard 
"C:\Program Files\R\R-2.4.1\bin\Rgui.exe" 
and paste the line into the box labeled "Type the location of the item"
4. Click the Next> button.
5. Type 'QRMBook' (without any quotation marks) into the box labeled "Type
a name of this shortcut".  Then click the Finish button.
6. Find the shortcut labeled "QRMBook" on your desktop,right-click the icon,
and choose 'Properties'.
7. The 'Start in' box says "C:\Program Files\R\R-2.4.1\bin".  Modify it to
read "C:\Program Files\R\R-2.4.1\Users\QRMBook" (be sure to include the 
quotation marks).  Then click OK.

Generally, you may now open the program by double-clicking the QRMBook shortcut on
your desktop.  However, you will then have to type the command
> library(QRMlib)
into the R-console each time you start QRMBook to load the the QRM library package.  
To avoid that step each time you open the shortcut, perform one of the following
two actions.  Once you have finished either of these actions, you will not need
to load the QRMlib into the workspace each time you start it.  
i) Copy the file .Rprofile from the 
	C:\Program Files\R\R-2.4.1\src\library\QRMlib\inst\Chap2
folder (substitute your version of R for R-2.4.1 if necessary) into the new
	C:\Program Files\R\R-2.4.1\Users\QRMBook
folder your created in step 1 above; OR
ii) Perform steps 1-7 in ..Creating an .Rprofile file... below

After you have finished installing an .Rprofile file into your QRMBook folder, you
can begin running scripts available for most of the chapters in the QRM book.
You will find scripts which explain topics in each QRM Book chapter. The folders 
are C:\Program Files\R\R-2.4.1\library\QRMlib\Chap2, ...\Chap3, etc.  You may open
these scripts by choosing File | Open Script within R, and then moving to the 
appropriate Chapter script for the QRM Book.

###Insuring data availability ###########
The scripts use data included with the installation.  The following data files are located at
	C:\Program Files\R\R-2.4.1\library\QRMlib\data
subfolder. If you examine that folder you may see the data files are compressed in a file
named Rdata.zip.  You may extract the data into the folder if you wish to see the names of
each separate data file by using WinZip or PKZip.  The data files include:
cac40.R, danish.R, DJ.R, dji.R, ftse100.R, FXGBP.RAW.R, hsi.R, nasdaq.R, nikkei.R, smi.R,
sp500.R, spdata.R, spdata.raw.R, and xdax.R.

Data files must be loaded into your workspace before they can be used.  The appropriate
command is
	data(filename)
where filename is the name of one of the data files WITHOUT its R extension.  Hence use
	data(sp500)
to load the data from the file sp500.R into the workspace.

When you exit the R program, you will be asked whether to save the current workspace 
environment.  If you choose to do so, the data files which you have opened via
data(filename) calls will be saved in the workspace so you don't need to execute any
subsequent data(filename) calls to get previously-loaded data.  If you do not save the
workspace, you must execute data(filename) each time you open the QRMBook workspace.


####Creating an .Rprofile file in "C:\Program Files\R\R-2.4.1\users\QRMBook"#######
#If you are using a later version of R like 2.4.0, replace 2.4.1 with 2.4.0 (or
#whatever version you are using).
1. Copy the next twelve lines of code (starting with # and ending with } into the clipboard.
#This file uses .First and .Last functions rather than the 
#  options(defaultPackages()) method to set the libraries to be loaded when
# a project starts.  
.First <- function()
{
   library(QRMlib)
}

.Last <- function()
{
   detach(package:QRMlib)
}

2. Open Notepad: left-click the Start button, choose Run and type notepad into
the box. We will try to save a file named ".Rprofile".  Note there will  be
no letters prior to the '.' and the type of file is an "Rprofile" type spelled
with a capital R followed by all small letters.
3. Paste the copied code into Notepad.
4. In Notepad, choose File | Save As from the menu.
5. In the resulting box, click the "Save as Type" drop-down box and choose "All Files".
(We are NOT saving as a .txt type.)
6. Paste the path "C:\Program Files\R\R-2.4.1\users\QRMBook\.Rprofile" into the 
File name box.
7. Click the Save button.


Work through the scripts in the Chap2 - Chap8 subdirectory.
These are organized according to chapters of the QRM book. 
There is no separate manual for the QRMlib library.  However, see the QRMlib.pdf
file in the Docs subfolder for a description of topics covered in the QRMlib.