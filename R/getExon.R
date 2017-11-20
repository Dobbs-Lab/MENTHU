#' getExon
#' 
#' Function that fetches exon data from a GenBank accesion ID
#' 
#' @param accession GenBank accession ID
#' @param percent Percentage of exons from 5' to 3' to fetch
#' @param wiggle If TRUE adds wigRoom nt of flanking intron to each exonSeq 
#' 
#' @result List of exon start and end site and DNA sequences
#' 
#' @examples 
#' 
#' @export 
#' 
getExon <- function(accession, percent = 100, wiggle = TRUE, wigRoom = 39) {
  
  require(genbankr)
  require(Biostrings)
  
  # Transform accession string to readable format and parse GenBank data
  gba <- GBAccession(accession) 
  gb <- readGenBank(gba)
  seq <- getSeq(gb)[[1]] # fetch gene sequence
  exonInfo <- data.frame(slot(exons(gb),"ranges")) #fetch exon information
  
  # Variable initialization
  set <- NULL
  exonSeq <- NULL
  
  # Calculation of number of exons to extract from percent parameter
  numExons <- floor(percent/100*length(exonInfo$start))
  
  if(wiggle == TRUE) {
    # Generation of DNAStringSet of exon DNA sequences
    for (i in 1:numExons) {
      exStart <- exonInfo$start[i] - wigRoom # Add intron wiggle room
      exEnd <- exonInfo$end[i] + wigRoom # Add intron wiggle room
      set <- c(set, seq[exStart:exEnd])
      exonSeq <- DNAStringSet(set)
    }
    temp <- exonInfo[1:numExons,]
    start <- temp$start - wigRoom
    end <- temp$end + wigRoom
    width <- temp$width + wigRoom*2
    table <- data.frame(start,end,width)
    return(list(table,exonSeq,seq))
  } else {
    # Generation of DNAStringSet of exon DNA sequences
    for (i in 1:numExons) {
      exStart <- exonInfo$start[i]
      exEnd <- exonInfo$end[i]
      set <- c(set, seq[exStart:exEnd])
      exonSeq <- DNAStringSet(set)
    }
    return(list(exonInfo[1:numExons,],exonSeq,seq))
  }
}
