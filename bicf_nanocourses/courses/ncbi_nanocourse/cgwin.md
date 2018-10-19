# Installation Instructions for CGWin

# Downloading Cygwin
* Download Cygwin from https://www.cygwin.com/  Select the setup-x86_64.exe file.
* Save the setup-x86_64.exe file to the Desktop on your laptop.

## Installing the Cygwin main program
* Double-click on the setup-x86_64.exe file on your Desktop. The Cygwin Net Release Setup Program will appear.  Select "Next."
* Select "Install from Internet", click "Next."
* Change the install directory from (C:\cygwin64) to a Cygwin directory on your Desktop (C:\Users\your_user_name\Desktop\Cygwin). If the folder does not exist, it will be created.  Select the "Just me" button is selected and click "Next".
* The dialog box asks for a local package directory.  Type Packages and click "Next". If the folder does not exist, it will be created.
* For the internet connection, leave "Direct Connection" selected and click "Next".
* For the Download site, select the first entry: Cygwin.mirror.constant.com selected and click "Next."
* Adding the Perl package, nano and unzip
* The Next window is for selecting packages. EDirect requires the complete Perl package and the unzip utility. The nano screen editor is helpful but not essential.  
* Next to the line for Perl, click on the word "default."  It should now display "install".  In the search bar at the top, enter nano.  
* The list will display Debug and Editors.  Click on Editors and an entry for nano will appear. If it says skip, click on the word skip and it should change to a version number (2.8.2-1)
* In the search bar, delete nano and type unzip. There will be two selections, Archive and Debug.
* Click on Archive and an entry for unzip will appear.  If it says "skip", click on "skip." A version number will appear (6.0-16).
* Now, click on "Next" in the lower right corner of the window.
* A message regarding resolving dependencies will appear.  Select "Next".

The installation should require approximately 15 minutes.

* Once installation is complete, a prompt for creating icons is displayed.  You'll probably want at least an icon on the Desktop. 
* Select "Finish."  Click on the Cygwin icon on the desktop or start menu.
* Now install [edirect](https://dataguide.nlm.nih.gov/edirect/install.html#edirect-installation)

Following all the directions, the entire process should require less than 30 minutes.
