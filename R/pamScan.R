#' pamScan
#' 
#' Function than searches a list of PAM sites in list of DNA targets
#' 
#' @param pamList list of PAMs and distance to cut site to screen for
#' @param targetList list of DNA targets to screen
#' @param exonStarts list of exon-start sites
#' 
#' @result list of PAM start sites for all PAMs in pamList in both 
#' strands for all targets queried
#' 
#' @examples 
#' 
#' @export 
#' 
pamScan <- function(pamList, targetList, exonStarts = NULL, findCut = FALSE, type = NULL) {
  
  require(Biostrings)
  
  # Variable definition and initialization
  numPAM <- length(pamList[[1]]) # Number of PAMs in pamList
  pamSites <- vector("list", numPAM)
  numTarget <- length(targetList) # Number of DNA targets to screen
  
  for (i in 1:numPAM) { # Iterate for each PAM
    
    # Variable initialization
    ppam <- DNAString(pamList[[1]][i]) # Transform char to DNAString
    npam <- reverseComplement(ppam) # Find PAM sequence in - DNA strand
    pPAM <- vector("list", numTarget)
    nPAM <- vector("list", numTarget)
    PAMs <- vector("list", 2)
    
    for (j in 1:length(targetList)) { # Iterate for each DNA target
      
      # Variable initialization
      target <- targetList[[j]]
      ppamSites <- NULL
      npamSites <- NULL
      
      # Find ppam and npam info in target and make df
      pdf <- data.frame(slot(matchPattern(ppam,target,fixed = FALSE),"ranges")) 
      ndf <- data.frame(slot(matchPattern(npam,target,fixed = FALSE),"ranges"))
      
      # Fetch start positions for ppam and npam sequences in DNA target 
      if (findCut == FALSE) {
        for (k in 1:length(pdf$start)) {
          ppamSites[k] <- pdf$start[k]
        }
        for (k in 1:length(ndf$start)) {
          npamSites[k] <- ndf$start[k]
        }
      } else {
        for (k in 1:length(pdf$start)) {
          if (pamList[[2]][i] < 0) {
            ppamSites[k] <- pdf$start[k] + pamList[[2]][i] - 1
          } else {
            ppamSites[k] <- pdf$start[k] + nchar(ppam) + pamList[[2]][i] - 1
          }
        }
        for (k in 1:length(ndf$start)) {
          if (pamList[[2]][i] < 0) {
            npamSites[k] <- ndf$start[k] + nchar(ppam) - pamList[[2]][i] - 1
          } else {
            npamSites[k] <- ndf$start[k] - pamList[[2]][i] - 1
          }
        }
      }
      
      # Add ppamStes and npamSites to memory
      pPAM[[j]] <- ppamSites
      nPAM[[j]] <- npamSites
      
      # Correct positions by localization in gene and not exon
      if (is.null(exonStarts) == FALSE) {
        pPAM[[j]] <- pPAM[[j]] + exonStarts[j] - 1
        nPAM[[j]] <- nPAM[[j]] + exonStarts[j] - 1
      }
    }
    
    # Name exons
    names(pPAM)[1:numTarget] <- paste0("Exon ", 1:(numTarget))
    names(nPAM)[1:numTarget] <- paste0("Exon ", 1:(numTarget))
    
    # Name PAMs by strand
    PAMs <-list(pPAM, nPAM)
    names(PAMs)[1:2] <- paste0(pamList[[1]][i], c(": Positive strand", ": Negative strand"))
    
    # Name exons by PAM
    pamSites[[i]] <- PAMs
    names(pamSites)[i] <- pamList[[1]][i]
  }
  return(pamSites)
}
