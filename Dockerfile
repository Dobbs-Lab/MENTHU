FROM openanalytics/r-base

MAINTAINER Carla Mann "genesculptsuitehelp@gmail.com"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0 \
    libxml2-dev

# install package dependencies for MENTHU
#RUN R -e "install.packages(c('shiny', 'XML', 'shinyjs', 'rhandsontable', 'plyr', 'stringr', 'stringi', 'rentrez', 'rlist', 'DT', 'devtools', 'httpuv', 'httr'), repos='https://cloud.r-project.org/')" -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("Biostrings")' -e 'devtools::install_github("rstudio/shiny-incubator")'
RUN R -e "install.packages(c('shiny', 'XML', 'xml2', 'shinyjs', 'rhandsontable', 'plyr', 'stringr', 'stringi', 'rentrez', 'rlist', 'DT', 'devtools', 'curl', 'plyr', 'jsonlite', 'httr'), repos='https://cloud.r-project.org/')" -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("Biostrings")' -e 'devtools::install_github("rstudio/shiny-incubator")'

# Copy MENTHU to image
RUN mkdir /root/menthu/
COPY / /root/menthu

COPY Rprofile.site /usr/lib/R/etc/

#Expose this port
EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/menthu/', port = 3838)"]
