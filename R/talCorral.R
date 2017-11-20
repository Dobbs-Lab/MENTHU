#' talCorral
#' 
#' Function than searches a list of PAM sites in list of DNA targets
#' 
#' @param talSeq T to A output sequence from talpal function
#' @param armin minimum allowable TALEN arm length
#' @param armax maximum allowable TALEN arm length
#' @param spamin minimum allowable spacer length
#' @param spamax maximum allowable spacer length
#' 
#' @result table of possible TALEN pairs
#' 
#' @examples 
#' 
#' @export 
#' 
talCorral <- function(talSeq, armin = 15, armax = 18, spamin = 14, spamax = 16) {
  
  numSeq <- length(talSeq)
  output <-vector("list", numSeq)
  
  # Remove 5'Ts from sequence
  for (i in 1:numSeq) {
    talSeq[[i]] <- substr(talSeq[[i]], 2, nchar(talSeq[[i]]) - 1)
  }
  
  for (k in 1:numSeq) {
    
    # Variable initialization
    tempL <- NULL
    tempR <- NULL
    spacer <- NULL
    talL <- NULL
    talR <- NULL
    spa <- NULL
    specs<- NULL
    
    # Find all possible TALEN arms
    for (i in armin:armax) {
      tempL[i-(armin - 1)] <- substr(talSeq[[k]], 1, i)
      tempR[i-(armin - 1)] <- substr(talSeq[[k]], nchar(talSeq[[k]]) - i + 1, nchar(talSeq[[k]]))
    }
    
    # Select for pairs that comply with spacer constraints
    for (i in 1:length(tempL)) {
      for (j in 1:length(tempR)) {
        spa <- substr(talSeq[[k]],(nchar(tempL[i]) + 1),(nchar(talSeq[[k]]) - nchar(tempR[j])))
        if (nchar(spa) >= spamin && nchar(spa) <= spamax) {
          spacer <- c(spacer,spa)
          talL <- c(talL,tempL[i])
          talR <- c(talR,tempR[j])
        }
      }
    }
    
    # Calculate arm and spacer lengths
    df <- cbind(nchar(talL), nchar(spacer), nchar(talR))
    for (i in 1: length(talL)) {
      specs <- c(specs, paste0(df[i,1], "/", df[i,2], "/", df[i,3]))
    }
    
    # Generate output table
    table <- NULL
    table <- cbind(talL, spacer, talR, specs)
    colnames(table) <- c("Left arm sequence", "Spacer sequence", "Right arm sequence", "Left/Spacer/Right length (nt)")
    output[[k]] <- table
  }
  
  return(output)
}
