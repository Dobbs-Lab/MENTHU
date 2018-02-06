# MENTHU
This is a repository for the MENTHU knockout site recommender.

You can run MENTHU online through a point-and-click website here: http://ll-g2f.gdcb.iastate.edu/menthu/

**If you already have R and/or RStudio installed, you can click [here](https://github.com/Dobbs-Lab/MENTHU/README.md#running-menthu-locally) to immmediately start running MENTHU locally.**

# How to Run MENTHU Locally
You will need to have the ability to install software on the computer you are using to run MENTHU locally; this may require administrator privileges. 

**If you are having issues running MENTHU locally, please check the 'Troubleshooting' section before requesting help.**


## Download and Install R
MENTHU requires the latest version of R in order to run offline. 

You can download R here: https://mirror.las.iastate.edu/CRAN/

Once you have downloaded the R installer, run it to install R. You may be required to enter administrator credentials; if you do not have these credentials, talk to your institution's IT department to have them install the software.

If you need additional help installing R, please check the installation instructions for your operating system:

Windows:    https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Windows

Mac OS:     https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-macOS

Linux/Unix: https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Unix_002dalikes


## Download and Install RStudio (optional)
MENTHU does not require the use of the RStudio development environment, but if you are interested in examining or modifying the MENTHU code, we recommend you do so in RStudio. 

You can download RStudio for free here: https://www.rstudio.com/products/rstudio/download/#download

After downloading the RStudio installer, follow the installation instructions. If you have both R and RStudio installed, you should only do the following steps in RStudio.


## Running MENTHU locally
You can run this RShiny web app in R (or RStudio) by opening up an R or RStudio session.

You can copy and paste the code blocks below into your R/RStudio console to run them.

(Please read [this discussion](https://www.lifehacker.com.au/2016/05/be-careful-when-you-copy-and-paste-code-from-the-internet/) on why you should not generally copy/paste/run code directly from internet pages to your console; you're probably safe on GitHub, but you should not get in the habit of copying and pasting code directly into your console or terminal windows. Paste what you've copied into NotePad or a similar program, and then copy/paste from there to your console once you know the code is safe.)


### Run this code ONLY THE FIRST TIME you run this tool on a computer, or when you need to update these packages:**

```
#Install packages required to run MENTHU; you can also run this code to update these packages
#Install CRAN packages
install.packages(c("shiny", "shinyjs", "Rcpp", "plyr", "stringr", "stringi", "shinyTable", "rentrez", "rlist", "DT", "xlsx", "devtools", "rhandsontable"))

#Install 'Biostrings' package from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

#Install 'ShinyIncubator' and 'ShinyTable' from GitHub
devtools::install_github("rstudio/shiny-incubator", force = TRUE)
devtools::install_github("trestletech/shinyTable",  force = TRUE)
```

### Run this code every time you want to use the tool, including the first time:**

```
#Load Shiny in the R/RStudio Environment
library(shiny)

#Retrieve, load, and run MENTHU from GitHub
runGitHub("MENTHU", "Dobbs-Lab")
```
You're all set!


# Troubleshooting
Please check the list of issues below to see if your issue is solved; if you can't find your issue or the fixes below don't work, please contact us at GeneSculptSuiteHelp@gmail.com for aid in troubleshooting. **Please be aware that we can only support the most up-to-date, unmodified versions of MENTHU's code.**

## Before troubleshooting or asking for help, please follow these simple steps that will resolve many (if not most) problems:
### 1. Check that you're running the most up-to-date version of MENTHU's code (currently v1.1.0) 
You can check this by looking at the upper right hand corner of the MENTHU user interface. This will generally only be an issue if you have downloaded the MENTHU code from GitHub and are using ```runApp()``` to run locally; if that is the case, download the updated code and see if your problem persists.


### 2. Update the required MENTHU packages 
You can update all required packages by re-running the code [here]().


### 3. Update your R or RStudio installation 
You can update R by following these directions:

On Windows:

Open R, and run the following code:

```
#Install and load the 'installr' package
install.packages("installr")
library(installr)

#This command will open a prompt to guide you through updating your R installation
updateR()
```

On Mac, Linux, or Unix:



## Potential issues:

### 1. Failure to load 'rJava'
Symptom: When trying to run the code, you get an error that contains:

```
Error : .onLoad failed in loadNamespace() for 'rJava'
```

Explanation: 'rJava' is a package that is required for downloading your results file.  'rJava' allows R/RStudio to communicate with Java on your computer to open a window to download and save your results file. When you get the ```loadNamespace() for 'rJava'``` error, this means that R/RStudio can't communicate with the Java on your system.

#### A1. Possible Cause: Java is not installed
Solution: Download and install Java.

Windows:
Visit https://www.java.com/en/download/help/windows_manual_download.xml and follow the Java installation instructions.

Mac OS:
Visit https://java.com/en/download/help/mac_install.xml and follow the Java installation instructions.

Linux/Unix:
Visit https://www.java.com/en/download/help/linux_install.xml and follow the Java installation instructions.


#### B1. Possible Cause: Java is out of date
Explanation: If you know you have Java installed on your system but are still getting the ```loadNamespace() for 'rJava'``` error, it is possible that an update to MENTHU and/or the packages it depends on requires an updated version of Java; also, Java installations can sometimes get corrupted or otherwise messed up, and this can cause R to not recognize Java on your system.

Solution: Download and install the latest version of Java for your operating system. Follow the directions in the links in **'B1. Possible Cause: Java is not installed'** to install updated Java.


#### C1. Possible Cause: R/RStudio can't find your Java installation
Explanation: If you know you have the latest version of Java installed on your system but you're still getting this error, it is possible that R/RStudio doesn't know where Java is located on your system. 


Solution 1: Add your Java installation location to your environment $PATH

Explanation

Solution 2:

Windows:

Mac OS:

1. Open a terminal window.

2. Run this code in the terminal:

```
sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
```

Linux/Unix:

Try running this code so that R can find your Java installation


#### D1. Possible Cause: There is a mismatch in your Java bit-version and R/RStudio bit-version
Explanation: R and RStudio can be downloaded in 32- and 64-bit versions. The differences between the 32- and 64-bit versions are not important to troubleshooting; what *is* important is that on some systems, Java will not communicate with R (and thus, 'rJava' will not load) if you have the 64-bit version of R and the 32-bit version of Java, or vice-versa. 

Solution: Download the bit-version of Java that matches your R/RStudio installation.

1. Identify your R bit-version:


2. Identify your Java bit-version:


3. If there is a mismatch, download the Java bit-version that matches your R bit-version:



