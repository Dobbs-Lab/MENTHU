##########Contains functions for dealing with target site identification

#' pamScan
#' 
#' Function than searches a list of PAM sites in list of DNA targets
#' 
#' @param pamList list of PAM sites to screen for
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
pamScan <- function(pamList, cutDistList, targetList, exonList, exonStarts = NULL, findCut = FALSE, type = NULL, wiggle = TRUE, wigRoom = 39) {
	
	require(Biostrings)
	
	# Variable definition and initialization
	numPAM <- length(pamList) # Number of PAMs in pamList
	pamSites <- vector("list", numPAM)
	numTarget <- length(targetList) # Number of DNA targets to screen
	
	for (i in 1:length(pamList)) { # Iterate for each PAM
		
		# Variable initialization
		ppam <- DNAString(pamList[[i]]) # Transform char to DNAString
		npam <- Biostrings::reverseComplement(ppam) # Find PAM sequence in - DNA strand
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
					if(cutDistList[i] < 0){
						ppamSites[k] <- pdf$start[k] + cutDistList[i] - 1
						
					} else {
						ppamSites[k] <- pdf$start[k] + nchar(ppam) + cutDistList[i] - 1
					}
					
				}
				for (k in 1:length(ndf$start)) {
					if (cutDistList[i] < 0) {
						npamSites[k] <- ndf$start[k] + nchar(ppam) - cutDistList[i] - 1
					} else {
						npamSites[k] <- ndf$start[k] - cutDistList[i] - 1
					}
				}
			}
			
			
			# Add ppamSites and npamSites to memory
			pPAM[[j]] <- ppamSites
			nPAM[[j]] <- npamSites
			
			# Correct positions by localization in gene and not exon
			if(!is.null(exonStarts)) {
				if(wiggle && (exonStarts[j] - wigRoom > 0)){
					pPAM[[j]] <- pPAM[[j]] + exonStarts[j] - wigRoom - 1
					nPAM[[j]] <- nPAM[[j]] + exonStarts[j] - wigRoom - 1
				} else if (wiggle && (exonStarts[j] - wigRoom < 0)){
					pPAM[[j]] <- pPAM[[j]] + exonStarts[j] - 1
					nPAM[[j]] <- nPAM[[j]] + exonStarts[j] - 1
				} else {
					pPAM[[j]] <- pPAM[[j]] + exonStarts[j]
					nPAM[[j]] <- nPAM[[j]] + exonStarts[j]
				}
			}
		}
		
		# Name exons if there is a list...
		#if()
		names(pPAM)[1:numTarget] <- paste0("Exon ", exonList)
		names(nPAM)[1:numTarget] <- paste0("Exon ", exonList)
		
		# Name PAMs by strand
		PAMs <-list(pPAM, nPAM)
		names(PAMs)[1:2] <- paste0(pamList[[i]], c(": Positive strand", ": Negative strand"))
		
		# Name exons by PAM
		pamSites[[i]] <- PAMs
		names(pamSites)[i] <- pamList[[i]]
	}

	return(pamSites)
}


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
talPal <- function(targetList, findcut = FALSE, wiggle = TRUE, wigRoom = 39, range = FALSE, 
									 armin = 15, armax = 18, spamin = 14, spamax = 16, exonList, exonStarts = NULL) {
	
	require(stringr)
	
	targetList <- as.character(targetList) # Transforms DNAString to char
	
	# Variable definition and initialization
	numSites <- length(targetList) # Number of sites to target
	talSites <- vector("list", numSites)
	talSeqs  <- vector("list", numSites)
	tempL <- NULL
	tempR <- NULL
	spacer <- NULL
	left <- NULL
	right <- NULL
	specs <- NULL
	
	# Calculates range of allowable spacer and arm lengths
	if (range == TRUE) {
		low     <- 2 * armin + spamin # Lower bound of distance between 5'Ts
		high    <- 2 * armax + spamax # Upper bound of distance between 5'Ts
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
			
			# Correct positions by localization in gene and not exon, taking into account "wiggle room"
			if (is.null(exonStarts) == FALSE) {
				if(wiggle && (exonStarts[i] - wigRoom > 0)){
					talSites[[i]] <- talSites[[i]] + exonStarts[i] - wigRoom - 1
					
				} else if (wiggle && (exonStarts[i] - wigRoom <= 0)) {
					talSites[[i]] <- talSites[[i]] + exonStarts[i] - 1
					
				} else {
					talSites[[i]] <- talSites[[i]] + exonStarts[i]
				}
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
				if(wiggle && (exonStarts[i] - wigRoom > 0)){
					talSites[[i]] <- talSites[[i]] + exonStarts[i] - wigRoom - 1
				} else if (wiggle && (exonStarts[i] - wigRoom <= 0)) {
					talSites[[i]] <- talSites[[i]] + exonStarts[i] - 1
				} else {
					talSites[[i]] <- talSites[[i]] + exonStarts[i]
				}
			}
		}
	}
	
	# Name exons
	names(talSites)[1:numSites] <- paste0("Exon ", exonList)
	names(talSeqs)[1:numSites] <- paste0("Exon ", exonList)
	
	return(list(talSites,talSeqs))
}

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
					talL   <- c(talL,tempL[i])
					talR   <- c(talR,tempR[j])
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
