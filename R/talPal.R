#' talPal
#' 
#' Function than searches a list of PAM sites in list of DNA targets
#' 
#' @param targetList list of DNA targets to screen
#' @param range if true, allows TALEN arm and spacer length flexibility
#' @param armin minimum allowable TALEN arm length
#' @param armax maximum allowable TALEN arm length
#' @param spamin minimum allowable spacer length
#' @param spamax maximum allowable spacer length
#' @param exonStarts list of exon-start sites
#' 
#' @result list of TALEN left arm 5'T location
#' 
#' @examples 
#' 
#' @export 
#' 

talPal <- function(targetList, findcut = FALSE, range = FALSE, armin = 15, armax = 18, spamin = 14, spamax = 16, exonStarts = NULL) {
  
  require(stringr)
  
  targetList <- as.character(targetList) # Transforms DNAString to char
  
  # Variable definition and initialization
  numSites <- length(targetList) # Number of sites to target
  talSites <- vector("list", numSites)

  talSeqs <- vector("list", numSites)
  tempL <- NULL
  tempR <- NULL
  spacer <- NULL
  left <- NULL
  right <- NULL
  specs <- NULL

  
  # Calculates range of allowable spacer and arm lengths
  if (range == TRUE) {
    low <- 2*armin + spamin # Lower bound of distance between 5'Ts
    high <- 2*armax + spamax # Upper bound of distance between 5'Ts
    pattern <- paste0("(?=(T[ACGT]{",low,",",high,"}A))") # Search regex
    
    # Find and extract left 5'T index and sequence from targetList
    for (i in 1:numSites) {
      target <- targetList[[i]]
      m <- gregexpr(pattern, target, perl = TRUE)
      tempn <- length(m[[1]])
      
      # Captures non-capturing look ahead assertion
      m <- lapply(m, function(i) {
        attr(i,"match.length") <- attr(i,"capture.length")
        i
      })
      
      talSeqs[[i]] <- regmatches(target,m) 
      
      if (findcut == TRUE) {
        # Estimates position of cut site
        talSites[[i]] <- unlist(as.list(m[[1]][1:tempn] + ceiling(attributes(m[[1]])$capture.length/2) - 1))
      } else{
        talSites[[i]] <- unlist(as.list(unlist(m[[1]][1:tempn])))
      }
      
      # Correct positions by localization in gene and not exon
      if (is.null(exonStarts) == FALSE) {
        talSites[[i]] <- talSites[[i]] + exonStarts[i] - 1

      }
    }
  } else {
    # Default TALEN search, 15b arm lengths and 15b spacer
    for (i in 1:numSites) {
      target <- targetList[[i]]

      m <- gregexpr("(?=(T[ACGT]{45}A))",target, perl = TRUE)
      tempn <- length(m[[1]])
      
      # Captures non-capturing look ahead assertion
      m <- lapply(m, function(i) {
        attr(i,"match.length") <- attr(i,"capture.length")
        i
      })
      
      talSeqs[[i]] <- regmatches(target,m)
      
      if (findcut == TRUE) {
        # Estimates position of cut site
        talSites[[i]] <- m[[1]][1:tempn] + ceiling(attributes(m[[1]])$capture.length/2) - 1 
      } else{
        talSites[[i]] <- m[[1]][1:tempn]
      }
      
      # Correct positions by localization in gene and not exon
      if (is.null(exonStarts) == FALSE) {
        talSites[[i]] <- talSites[[i]] + exonStarts[i]
      }
    }
  }
  
  # Name exons
  names(talSites)[1:numSites] <- paste0("Exon ", 1:(numSites))

  names(talSeqs)[1:numSites] <- paste0("Exon ", 1:(numSites))
  
  return(list(talSites,talSeqs))
}
