# MENTHU
# A repository for MENTHU knockout site recommender

#You can run this shiny web app in RStudio (or R) by opening up an R or RStudio session, and running the following code:

#Run first time you run this tool on a computer:

install.packages(c("shiny", "shinyjs", "Rcpp", "plyr", "stringr", "stringi", "shinyTable", "rentrez", "rlist", "DT", "xlsx", "devtools", "rhandsontable"))

source("https://bioconductor.org/biocLite.R")

biocLite("Biostrings")

devtools::install_github("rstudio/shiny-incubator")


#Run this code whenever you want to use the tool, including the first time:

library(shiny)

runGitHub("MENTHU", "Dobbs-Lab")


#You're all set!
