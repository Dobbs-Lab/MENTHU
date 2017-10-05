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
talPal <- function(targetList, range = FALSE, armin = 14, armax = 18, spamin = 14, spamax = 16, exonStarts = NULL) {
  
  require(stringr)
  
  targetList <- as.character(targetList) # Transforms DNAString to char
  
  # Variable definition and initialization
  numSites <- length(targetList) # Number of sites to target
  talSites <- vector("list", numSites)
  
  # Calculates range of allowable spacer and arm lengths
  if (range == TRUE) {
    low <- 2*armin + spamin # Lower bound of distance between 5'Ts
    high <- 2*armax + spamax # Upper bound of distance between 5'Ts
    pattern <- paste0("(?=T[ACGT]{",low,",",high,"}A)") # Search regex
    
    # Find pattern in targetList and extract left 5'T positions
    for (i in 1:numSites) {
      target <- targetList[[i]]
      templ <- gregexpr(pattern, target, perl = TRUE)
      tempn <- length(templ[[1]])
      talSites[[i]] <- templ[[1]][1:tempn]
      
      # Correct positions by localization in gene and not exon
      if (is.null(exonStarts) == FALSE) {
        talSites[[i]] <- talSites[[i]] + exonStarts[i]
      }
    }
  } else {
    # Default TALEN search, 15b arm lengths and 15b spacer
    for (i in 1:numSites) {
      target <- targetList[[i]]
      templ <- gregexpr("(?=T[ACGT]{45}A)",target, perl = TRUE)
      tempn <- length(templ[[1]])
      talSites[[i]] <- templ[[1]][1:tempn]
      
      # Correct positions by localization in gene and not exon
      if (is.null(exonStarts) == FALSE) {
        talSites[[i]] <- talSites[[i]] + exonStarts[i]
      }
    }
  }
  
  # Name exons
  names(talSites)[1:numSites] <- paste0("Exon ", 1:(numSites))
  
  return(talSites)
}
