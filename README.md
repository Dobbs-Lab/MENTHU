# MENTHU
This is a repository for the MENTHU knockout site recommender.

You can run MENTHU online through a web interface here: http://genesculpt.org/menthu/
 

### If you already have R and/or RStudio installed, you can jump to [here](https://github.com/Dobbs-Lab/MENTHU#3-run-menthu-locally) to immmediately start running MENTHU locally.

### If you are having issues running MENTHU locally, please check the [Troubleshooting](https://github.com/Dobbs-Lab/MENTHU#troubleshooting) section before requesting help.
 

# How to Run MENTHU Locally
You will need to have the ability to install software on the computer you are using to run MENTHU locally; this may require administrator privileges. 

[1. Download and Install R](https://github.com/Dobbs-Lab/MENTHU#1-download-and-install-r)

[2. Download and Install RStudio](https://github.com/Dobbs-Lab/MENTHU#2-download-and-install-rstudio-optional) (optional)

[3. Run MENTHU locally](https://github.com/Dobbs-Lab/MENTHU#3-run-menthu-locally)

[Troubleshooting](https://github.com/Dobbs-Lab/MENTHU#troubleshooting)

## [1. Download and Install R](#1-download-and-install-r)
MENTHU requires the latest version of R in order to run offline. 

Download R for your appropriate operating system:

Windows: https://mirror.las.iastate.edu/CRAN/bin/windows/

 You should select the "base" option, or click "install R for the first time".
 

Mac OS: https://mirror.las.iastate.edu/CRAN/bin/macosx/

 Scroll down to the "Files" section, and find the R pkg file that lists your operating system (El Capitan, Mavericks, Snow Leopard, etc). Select the R-3.x.x.pkg file corresponding to your system - pay special attention to the "Important" section under R-3.4.3.pkg if you have "El Capitan"; you may want to consider using R-3.3.3.pkg if you don't want to install additional tools to support running R 3.4.3 on "El Capitan".


Linux/Unix: https://mirror.las.iastate.edu/CRAN/bin/linux/

 Find your Unix distro, and folow the instructions in the directory.
 

Once you have downloaded the R installer, run it to install R. You may be required to enter administrator credentials; if you do not have these credentials, talk to your institution's IT department to have them install the software.


If you need additional help installing R, please check the installation instructions for your operating system:

Windows:    https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Windows

Mac OS:     https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-macOS

Linux/Unix: https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Unix_002dalikes



## [2. Download and Install RStudio](#2-download-and-install-rstudio-optional) (optional)
MENTHU does not require the use of the RStudio development environment, but if you are interested in examining or modifying the MENTHU code, we recommend you do so in RStudio. 

You can download RStudio for free here: https://www.rstudio.com/products/rstudio/download/#download

After downloading the RStudio installer, follow the installation instructions. If you have both R and RStudio installed, you should only do the following steps in RStudio.



## [3. Run MENTHU locally](#run-menthu-locally)
You can run this RShiny web app in R (or RStudio) by opening up an R or RStudio session.

You can copy and paste the code blocks below into your R/RStudio console to run them.

(Please read [this discussion](https://www.lifehacker.com.au/2016/05/be-careful-when-you-copy-and-paste-code-from-the-internet/) on why you should not generally copy/paste/run code directly from internet pages to your console; you're probably safe on GitHub, but you should not get in the habit of copying and pasting code directly into your console or terminal windows. Paste what you've copied into NotePad or a similar program, and then copy/paste from there to your console once you know the code is safe.)



### [Run this code ONLY THE FIRST TIME you run this tool on a computer, or when you need to update these packages:](#install-code)

```
#Install packages required to run MENTHU; you can also run this code to update these packages

#Install CRAN packages
install.packages(c("shiny", "shinyjs", "Rcpp", "plyr", "stringr", "stringi", "shinyTable", 
                   "rentrez", "rlist", "DT", "xlsx", "devtools", "rhandsontable", "httr", "jsonlite", "xml2"))

#Install 'Biostrings' package from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

#Install 'ShinyIncubator' from GitHub
devtools::install_github("rstudio/shiny-incubator", force = TRUE)
```

### Run this code every time you want to use the tool, including the first time:

```
#Retrieve, load, and run MENTHU from GitHub
shiny::runGitHub("MENTHU", "Dobbs-Lab")
```

You're all set!



# [Troubleshooting](#troubleshooting)
Please check the list of issues below to see if your issue is solved; if you can't find your issue or the fixes below don't work, please contact us at GeneSculptSuiteHelp@gmail.com for aid in troubleshooting. 

**Please be aware that we can only support up-to-date and unmodified versions of MENTHU's code.**

## Solutions:
[Simple Fixes That Resolve Most Problems](https://github.com/Dobbs-Lab/MENTHU#simple-fix)

[Known (And Solved) Issues](https://github.com/Dobbs-Lab/MENTHU#known-and-solved-issues)

## [Before troubleshooting or asking for help, please follow these simple steps that will resolve many (if not most) problems:](#simple-fix)
[1. Check that you're running up-to-date MENTHU](https://github.com/Dobbs-Lab/MENTHU#1-up-to-date)

[2. Update Required MENTHU packages](https://github.com/Dobbs-Lab/MENTHU#2-update-the-required-menthu-packages)

[3. Update R and/or RStudio installations](https://github.com/Dobbs-Lab/MENTHU#3-update-your-r-or-rstudio-installation)


### [1. Check that you're running the most up-to-date version of MENTHU's code](#1-up-to-date)
The most current MENTHU version is v1.1.0.

You can check which version you are using by looking at the upper right hand corner of the MENTHU user interface. This will generally only be an issue if you have downloaded the MENTHU code from GitHub and are using ```runApp()``` to run locally; if that is the case, download the updated code and see if your problem persists.


### [2. Update the required MENTHU packages](#2-update-the-required-menthu-packages)
You can update all required packages by re-running the code [here](https://github.com/Dobbs-Lab/MENTHU#install-code).


### [3. Update your R or RStudio installation](#3-update-your-r-or-rstudio-installation)
You can check if you have the latest version of R by running the command ```R.Version()``` in R or RStudio.

Find the ```$version.string``` output, which will say "R version x.x.x (20yy-mm-dd)". The numbers replacing "x.x.x" are your R version; you should be running R 3.4.3 if you are on Windows or Linux/Unix, and either R 3.4.3 or R 3.3.3 for Mac OS X 10.9 (Mavericks) and above.


If your R is out of date, you can update your installation by following these directions:


On Windows:

Open R, and run the following code:

```
#Install and load the 'installr' package
install.packages("installr")
library(installr)

#This command will open a prompt to guide you through updating your R installation
updateR()
```

You can also follow the instructions under [1. Download and Install R](https://github.com/Dobbs-Lab/MENTHU#1-download-and-install-r), but the ```installr``` package is probably quicker and easier.


On Mac, Linux, or Unix:
Follow the instructions under [1. Download and Install R](https://github.com/Dobbs-Lab/MENTHU#1-download-and-install-r) for your operating system.



## Known (and Solved) Issues:
[1. Failure to load 'rJava'](https://github.com/Dobbs-Lab/MENTHU#1-failure-to-load-rjava)

### 1. Failure to load 'rJava'
Symptom: When trying to run the code, you get an error that contains:

```
Error : .onLoad failed in loadNamespace() for 'rJava'
```

Explanation: 'rJava' is a package that is required for downloading your results file.  'rJava' allows R/RStudio to communicate with Java on your computer to open a window to download and save your results file. When you get the ```loadNamespace() for 'rJava'``` error, this means that R/RStudio can't communicate with the Java on your system.

Solutions:

[A1. Possible Cause: Java is not installed](https://github.com/Dobbs-Lab/MENTHU#a1-possible-cause-java-is-not-installed)

[B1. Possible Cause: Java is out of date](https://github.com/Dobbs-Lab/MENTHU#b1-possible-cause-java-is-out-of-date)

[C1. Possible Cause: R/RStudio can't find your Java installation](https://github.com/Dobbs-Lab/MENTHU#c1-possible-cause-rrstudio-cant-find-your-java-installation)

[D1. Possible Cause: There is a mismatch in your Java bit-version and R/RStudio bit-version](https://github.com/Dobbs-Lab/MENTHU#d1-possible-cause-there-is-a-mismatch-in-your-java-bit-version-and-rrstudio-bit-version)

#### [A1. Possible Cause: Java is not installed](#a1-possible-cause-java-is-not-installed)
You can check if Java is correctly installed on your system by running ```java -version``` in a terminal window (this works on Windows, Mac OS, and Linux/Unix systems.) If you get an output, you have Java correctly installed. If you get an error such as "java: Command not found", you either do not have Java installed or the terminal can't find your Java installation.

Solution: Download and install Java.

Windows:
Visit https://www.java.com/en/download/help/windows_manual_download.xml and follow the Java installation instructions.

Mac OS:
Visit https://java.com/en/download/help/mac_install.xml and follow the Java installation instructions.

Linux/Unix:
Visit https://www.java.com/en/download/help/linux_install.xml and follow the Java installation instructions.


#### [B1. Possible Cause: Java is out of date](#b1-possible-cause-java-is-out-of-date)
Explanation: If you know you have Java installed on your system but are still getting the ```loadNamespace() for 'rJava'``` error, it is possible that an update to MENTHU and/or the packages it depends on requires an updated version of Java; also, Java installations can sometimes get corrupted or otherwise messed up, and this can cause R to not recognize Java on your system.

Solution: Download and install the latest version of Java for your operating system. Follow the directions in the links in [A1. Possible Cause: Java is not installed](https://github.com/Dobbs-Lab/MENTHU#a1-possible-cause-java-is-not-installed) to install updated Java.


#### [C1. Possible Cause: R/RStudio can't find your Java installation](#c1-possible-cause-rrstudio-cant-find-your-java-installation)
Explanation: If you know you have the latest version of Java installed on your system but you're still getting this error, it is possible that R/RStudio doesn't know where Java is located on your system. 


Solution: Add your Java installation location to your $PATH variable

Explanation: The $PATH variable is a collection of directions to locations on your computer where executable programs are stored. When a program that requires Java to run goes to use Java, it looks in the $PATH variable to find where Java is on your computer. If Java's location is not in the $PATH variable, then the program can't find it, and so thinks you don't have Java installed.

Follow the instructions here https://www.java.com/en/download/help/path.xml to add Java to your system $PATH.


#### [D1. Possible Cause: There is a mismatch in your Java bit-version and R/RStudio bit-version](#d1-possible-cause-there-is-a-mismatch-in-your-java-bit-version-and-rrstudio-bit-version)
Explanation: R and RStudio can be downloaded in 32- and 64-bit versions. The differences between the 32- and 64-bit versions are not important to troubleshooting; what *is* important is that on some systems, Java will not communicate with R (and thus, 'rJava' will not load) if you have the 64-bit version of R and the 32-bit version of Java, or vice-versa. 

Solution: Download the bit-version of Java that matches your R/RStudio installation.

1. Identify your R bit-version:

      In R or RStudio, type ```R.Version()``` in the console.

      There will be several lines of output; look for the one starting ```$arch```.

      For 32-bit versions, this value will be "i386" or "x86".

      For 64-bit versions, this value will be "x86_64" or "x64"


2. Identify your Java bit-version:

      Open a terminal window, and type ```java -version```.

      If you see something like

       java version "1.x.x_xxx"
       Java(TM) SE Runtime Environment <build 1.x.x_xxx-xxx>
       Java Hostpot(TM) Client VM (build xx.xxx-xxxx, mixed mode)
       
      you have 32-bit Java installed.


      If you see something like

       java version "1.x.x_xxx"
       Java(TM) SE Runtime Environment <build 1.x.x_xxx-xxx>
       Java HotSpot<TM> 64-Bit Server VM <build xx.xxx-xxx, mixed mode>

      then you have 64-bit Java installed.


3. If there is a mismatch, download the Java bit-version that matches your R bit-version by following the directions under [A1. Possible Cause: Java is not installed](https://github.com/Dobbs-Lab/MENTHU#a1-possible-cause-java-is-not-installed)
.

