#' getExon
#' 
#' Function that extracts a subset of a gene sequence centered at position
#' 
#' @param sequence genetic sequence input
#' @param position position to center window at
#' @param winSize set sequence length to return
#' 
#' @result List of exon start and end site and DNA sequences
#' 
#' @examples 
#' 
#' @export 
#' 
window <- function(sequence, position, winSize = 80) {

return(sequence[(position-((winSize + winSize%%2)/2 - 1)):(position+((winSize + winSize%%2)/2))])
  
}
