#' calculateMenthu2
#'
#' @param inData  A DNA sequence window/context
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight  A weight factor for caculating pattern score; default is 20, as in Bae paper
#' @param maxdbm  Maximum distance between microhomologies
#'
#' @return

#For each sequence in the list, calculate the slope competition
calculateMenthu2 <- function(inData, cutSite = -1, weight = 20, maxdbm = 5){
	
	# If the cut site isn't specified, assume it is in the middle of the input
	if(cutSite == -1){
		cutSite <- floor(nchar(inData) / 2)
	}
	
	# Get ordered pattern scores for every deletion pattern in the wildtype sequence
	patternScoreDF  <- baeRevised(inData, cutSite, weight)
	
	if(nrow(patternScoreDF) > 0){
		#Sort frame by pattern score
		patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore), ]
		
		# Filter out MHs shorter than 3 from patternScoreDF
		threePlus      <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3),]
		
		# Filter out entries with top predicted deletion with more than 5bp separation between MHs
		if(threePlus$delLength[1] - nchar(patternScoreDF$microhomology)[1] <= maxdbm) {
			critOne  <- TRUE
		} else {
			critOne  <- FALSE
		}
		
		# Generate output values
		# If no entries comply with criterion (no 3+bp mh w/in 5bp), set ratio to zero and fShift to NA
		if(critOne == FALSE) { 
			ratio    <- 0
			fShift   <- "NA"
			
			# If only one entry complies with criterion, set an infinite ratio and determine fShift
		} else if(nrow(threePlus) == 1) { 
			ratio    <- Inf
			fShift   <- (if(threePlus$delLength[1] %% 3 == 0){"No"} else {"Yes"})
			
		} else { 
			# Calculation of MENTHU score 2.0	
			ratio    <- threePlus$patternScore[1] / threePlus$patternScore[2]
			
			#To prevent choking on 0/0 instances (no 3+bp microhomologies in anything), set ratio to zero and fShift to NA
			if(is.nan(ratio)){ 
				ratio  <- 0
				fShift <- "NA"
				
			} else {
				fShift <- (if(threePlus$delLength[1] %% 3 == 0){"No"} else {"Yes"})
			}
			
		}
		
		# Organize output values as data frama
		outFrame <- data.frame(seq              = inData,
													 menthuScore      = ratio,
													 frameShift       = fShift,
													 topDel           = threePlus$delSeq[1],
													 topMH            = threePlus$microhomology[1],
													 stringsAsFactors = FALSE)
	} else {
		#If there are no microhomologies, return an empty-ish frame
		outFrame <- data.frame(seq              = "", 
													 menthuScore      = -1, 
													 frameShift       = "NA", 
													 topDel           = "",
													 topMH            = "",
													 stringsAsFactors = FALSE)
	}
	
	return(outFrame)
}


#' calculateSlopeCompetition
#'
#' @param inData  A DNA sequence/window/context
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight  A weight factor for caculating pattern score; default is 20, as in Bae paper
#' @param top     The number of entries to consider when generating the linear regression
#'
#' @return

calculateSlopeCompetition <- function(inData, cutSite = -1, weight = 20, top = 10){
	#If the cut site isn't specified, assume it is in the middle of the input
	if(cutSite == -1){
		cutSite <- floor(nchar(inData) / 2)
	}
	
	#Get the pattern scores for every deletion pattern in the wildtype sequence
	patternScoreDF <- baeRevised(inData, cutSite, weight, mh = 2)
	patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore), ]
	
	if(nrow(patternScoreDF) > 0){
		#Determine if there is a frameshift for the top-scoring deletion pattern
		fShift <- (if(patternScoreDF[1, 3] %% 3 == 0){"No"} else {"Yes"})
		
		#Subset the out of frame scores
		outOfFrameInst <- patternScoreDF[which(patternScoreDF$delLength %% 3 != 0), ]
		
		#Calculate the microhomology score according to Bae algorithm
		mhScore <- sum((nchar(patternScoreDF$microhomology) +
											stringr::str_count(toupper(patternScoreDF$microhomology), "G") +
											stringr::str_count(toupper(patternScoreDF$microhomology), "C")) *
									 	(1 / exp((patternScoreDF$delLength) / weight)) * 100)
		
		#Calculate the out-of-frame score according to Bae algorithm
		outOfFrameScore <- (sum((nchar(outOfFrameInst$microhomology) +
														 	stringr::str_count(toupper(outOfFrameInst$microhomology), "G") +
														 	stringr::str_count(toupper(outOfFrameInst$microhomology), "C")) *
															(1 / exp((outOfFrameInst$delLength) / weight)) * 100) / mhScore) * 100
		
		
		###This subsets to only 3bp or greater microhomologies#########
		threePlus <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3), ]
		
		#If there are at least 10 3bp microhomologies, take the top 10; otherwise, use as many as are available
		if(nrow(threePlus) < top){
			end3 <- nrow(threePlus)
			
		} else {
			end3 <- top
			
		}
		
		#Generate list of top 10 3bp+ microhomologies
		top103Plus <- threePlus[1:end3, ]
		
		#Perform linear regression on the top 10
		linModel3 <- lm(top103Plus$patternScore ~ seq(1:end3), data = top103Plus)
		
		#Create data frame to hold information
		tempFrame <- data.frame(seq                = inData,
														microhomologyScore = mhScore,
														outOfFrameScore    = outOfFrameScore,
														slopeMH3Plus       = linModel3$coefficients[2],
														frameShift         = fShift,
														topDel             = top103Plus$seq[1],
														stringsAsFactors   = FALSE)
		
		
	} else {
		#Create empty data frame if there are no pattern scores
		tempFrame <- data.frame(seq                = "",
														microhomologyScore = 0,
														slopeMH3Plus       = 0,
														frameShift         = "NA",
														stringsAsFactors   = FALSE)
	}
	
	rownames(tempFrame) <- c()
	return(tempFrame)
}


#' baeRevised
#' 
#' @param seq         DNA sequence
#' @param cutPosition Position to be cut in the DNA sequence
#' @param weight      deletion length weight factor
#'
#' @return
#' @export
#'
#' @examples

baeRevised <- function(seq, cutPosition, weight = 20.0, mhL = 3){
	#Split the sequence into upstream and downstream of the cut site
	upstream   <- substring(seq, 1,                    nchar(seq) / 2)
	downstream <- substring(seq, (nchar(seq) / 2) + 1, nchar(seq))
	
	#Get all common substrings between upstream and downstream sections of sequence
	commonSubstrings <- allsubstr(upstream, downstream, mhL)
	
	#Check if there is a 3bp or greater MH
	if(length(commonSubstrings) > 0){ 
		
		#Find all upstream matches to common substrings
		matchLocsUp       <- stack(sapply(commonSubstrings, Biostrings::gregexpr2, text = upstream))
		matchLocsUp$ind   <- as.character(matchLocsUp$ind)
		
		#Find all downstream matches to common substrings
		matchLocsDown     <- stack(sapply(commonSubstrings, Biostrings::gregexpr2, text = downstream))
		matchLocsDown$ind <- as.character(matchLocsDown$ind)
		
		#Find all the starting and ending locations of the common substrings in the upstream portion
		upStart   <- matchLocsUp$values
		upEnd     <- matchLocsUp$values + nchar(matchLocsUp$ind) - 1
		
		#Do the same for the downstream
		downStart <- matchLocsDown$values + nchar(upstream)
		downEnd   <- matchLocsDown$values + nchar(matchLocsDown$ind) + nchar(upstream)
		
		#Store microhomologies in data frames
		mhDFup <- data.frame(seq              = matchLocsUp$ind,
												 upStart          = upStart,
												 upEnd            = upEnd,
												 stringsAsFactors = FALSE
		)
		
		mhDFdown <- data.frame(seq              = matchLocsDown$ind,
													 downStart        = downStart,
													 downEnd          = downEnd,
													 stringsAsFactors = FALSE
		)
		
		#Order the data frames
		mhDFup   <-   mhDFup[order(nchar(mhDFup$seq),   mhDFup$seq),   ]
		mhDFdown <- mhDFdown[order(nchar(mhDFdown$seq), mhDFdown$seq), ]
		
		#Merge the two data frames
		mhDF <- suppressMessages(plyr::join(mhDFup, mhDFdown))
		
		#Get the lengths of all the deletions
		leng <- mhDF$downEnd - mhDF$upEnd - 1
		
		#Get the microhomologies
		mh   <- mhDF$seq
		
		#Create the deletion sequences
		delSeqs <- unlist(lapply(1:nrow(mhDF), function(x) paste0(substring(seq, 1, mhDF$upEnd[x]),
																															paste(rep('-', leng[x]), collapse = ''),
																															substring(seq, mhDF$downEnd[x], nchar(seq)))))
		
		# Create the sequence that would be observed from the deletion
		delPattern <- gsub('-', '', delSeqs, fixed = TRUE)
		
		
		oof     <- unlist(lapply(leng, function(x) (if(x %% 3 != 0){1} else {0})))
		
		# Create a data frame to hold the information
		delFrame <- data.frame(seq               = rep(seq, length(delSeqs)),
													 deletedSeqContext = delSeqs,
													 delPattern        = delPattern,
													 microhomology     = mh, 
													 startDel          = mhDF$upEnd + 1, 
													 endDel            = mhDF$downEnd, 
													 mhStart1          = mhDF$upStart, 
													 mhEnd1            = mhDF$upEnd,
													 mhStart2          = mhDF$downStart,
													 mhEnd2            = mhDF$downEnd,
													 deletedSeq        = delSeqs,
													 delLength         = leng, 
													 mhLength          = nchar(mhDF$seq),
													 GC                = (stringr::str_count(mh, 'G') + stringr::str_count(mh, 'C')) / nchar(mh),
													 outOfFrame        = oof,
													 patternScore      = (100 * (round(1 / exp((leng) / weight), 3)) * 
													 										 	(stringr::str_count(mh, 'G') + stringr::str_count(mh, 'C') + nchar(mh))),
													 stringsAsFactors  = FALSE)
		
		# Order the data frame by microhomology length (descending), and then by pattern score (descending)
		delFrameOrd <- delFrame[order(-nchar(delFrame$microhomology), -delFrame$patternScore),]
		
		# Identify if a row in the data frame produces a deletion pattern produced by a longer microhomology
		# This eliminates redundant microhomologies which are actually part of a longer, better-scoring microhomology
		dupes <- sapply(1:nrow(delFrameOrd), function(x) unlist(sapply(1:x, function(y){
			if(x != y){
				if(delFrameOrd$delPattern[x] == delFrameOrd$delPattern[y]){
					TRUE
				} else {
					FALSE
				}
			} else {
				FALSE
			}
		})))
		
		# Drop the duplicates
		dupeDrop     <- sapply(dupes, function(x) any(x))
		dupeDropList <- which(dupeDrop)
		delFrameOrdDeDupe <- delFrameOrd[-dupeDropList,]
		delFrame          <- delFrameOrdDeDupe[order(-delFrameOrdDeDupe$patternScore), ]
		
		#Create the deletion/patternscore data frame
		psDF <- data.frame(seq              = delFrame$seq,
											 microhomology    = delFrame$microhomology,
											 delLength        = delFrame$delLength,
											 patternScore     = delFrame$patternScore,
											 delSeq           = delFrame$deletedSeqContext,
											 stringsAsFactors = FALSE
		)
		
	} else {
		
		#Create empty return frame if no MHs detected
		psDF <- data.frame(seq              = as.character(),
											 microhomology    = as.character(),
											 delLength        = as.numeric(),
											 patternScore     = as.numeric(),
											 delSeqs          = as.character(),
											 stringsAsFactors = FALSE)
	}
	
	return(psDF)
}


#' allsubstr
#'
#' @param upstream   upstream section of DNA sequence 
#' @param downstream downstream section of DNA sequence
#'
#' @return
#' @export
#'
#' @examples

allsubstr <- function(upstream, downstream, mh = 3){
	#Create empty string to hold upstream strings
	upstreamStrings   <- list()
	#Create empty string to hold downstream strings
	downstreamStrings <- list()
	
	#Generate all possible substrings of length >= 3 in the upstream and downstream strings
	for(i in mh:nchar(upstream)){
		upstreamStrings   <- c(upstreamStrings,   unique(substring(upstream,   1:(nchar(upstream)   - i + 1), i:nchar(upstream))))
		downstreamStrings <- c(downstreamStrings, unique(substring(downstream, 1:(nchar(downstream) - i + 1), i:nchar(downstream))))
	}
	
	#Find the intersection (common strings) between the upstream and downstream section
	commonStrings <- intersect(unlist(upstreamStrings), unlist(downstreamStrings))
	
	return(commonStrings)
} 
