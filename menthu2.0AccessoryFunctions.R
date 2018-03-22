#' calculateMenthu2
#'
#' @param inData A DNA sequence/window/context
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight A weight factor for caculating pattern score; default is 20, as in Bae paper
#' @param mindbm Maximum distance between microhomologies
#'
#' @return


#For each sequence in the list, calculate the slope competition
calculateMenthu2 <- function(inData, cutSite = -1, weight = 20, maxdbm = 5){
	
	# If the cut site isn't specified, assume it is in the middle of the input
	if(cutSite == -1){
		cutSite <- floor(nchar(inData)/2)
	}
	
	# Get ordered pattern scores for every deletion pattern in the wildtype sequence
	patternScoreDF <- calculateBae(inData, cutSite, weight)
	patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore),]
	
	if(nrow(patternScoreDF) > 0){
		
		
		# Filter out MHs shorter than 3 from patternScoreDF
		threePlus <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3),]
		
		# Filter out entries with top predicted deletion with more than 5bp separation between MHs
		if (threePlus$delLength[1] - nchar(patternScoreDF$microhomology)[1] <= maxdbm) {
			critOne <- TRUE
		} else {
			critOne <- FALSE
		}
		
		# Generate output values
		if (critOne == FALSE) { # If no entries comply with criterion (no 3+bp mh w/in 5bp)
			ratio  <- 0
			fShift <- "NA"
		} else if (nrow(threePlus) == 1) { # If only one entry complies with criterion
			ratio  <- Inf
			fShift <- (if(threePlus[1, 3] %% 3 == 0){"No"} else {"Yes"})
		} else { # Calculation of MENTHU score 2.0
			ratio <- threePlus[1,4]/threePlus[2,4]
			
			if(is.nan(ratio)){ #To prevent choking on 0/0 instances (no 3+bp microhomologies in anything)
				ratio  <- 0
				fShift <- "NA"
			} else {
				fShift <- (if(threePlus[1, 3] %% 3 == 0){"No"} else {"Yes"})
			}
			
		}
		
		# Organize output values as data frama
		outFrame <- data.frame(seq         = inData,
													 menthuScore = ratio,
													 frameShift  = fShift,
													 stringsAsFactors = FALSE)
	} else {
		outFrame <- data.frame(seq = "", menthuScore = -1, frameShift = "NA", stringsAsFactors = FALSE)
	}
	
	return(outFrame)
}

#' calculateSlopeCompetition
#'
#' @param inData A DNA sequence/window/context
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight A weight factor for caculating pattern score; default is 20, as in Bae paper
#' @param top The number of entries to consider when generating the linear regression
#'
#' @return


#For each sequence in the list, calculate the slope competition
calculateSlopeCompetition <- function(inData, cutSite = -1, weight = 20, top = 10){
	lengthWeight <- weight
	
	#for(a in 1:nrow(inData)){
	if(cutSite == -1){
		#If the cut site isn't specified, assume it is in the middle of the input
		cutSite <- floor(nchar(inData)/2)
	}
	
	#Get the pattern scores for every deletion pattern in the wildtype sequence
	patternScoreDF <- calculateBae(inData, cutSite, lengthWeight)
	patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore),]
	
	if(nrow(patternScoreDF) > 0){
		fShift <- (if(patternScoreDF[1, 3] %% 3 == 0){"No"} else {"Yes"})
		
		#Subset the out of frame scores
		outOfFrameInst <- patternScoreDF[which(patternScoreDF$delLength %% 3 != 0),]
		
		#Calculate the microhomology score according to Bae algorithm
		mhScore <- sum((nchar(patternScoreDF$microhomology) +
											stringr:::str_count(toupper(patternScoreDF$microhomology), "G") +
											stringr:::str_count(toupper(patternScoreDF$microhomology), "C")) *
									 	(1/exp((patternScoreDF$delLength)/lengthWeight)) * 100)
		
		#Calculate the out-of-frame score according to Bae algorithm
		outOfFrameScore <- (sum((nchar(outOfFrameInst$microhomology) +
														 	stringr:::str_count(toupper(outOfFrameInst$microhomology), "G") +
														 	stringr:::str_count(toupper(outOfFrameInst$microhomology), "C")) *
															(1/exp((outOfFrameInst$delLength)/lengthWeight)) * 100) / mhScore) * 100
		
		
		###This subsets to only 3bp or greater microhomologies#########
		threePlus <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3),]
		
		#If there are at least 10 3bp microhomologies, take the top 10; otherwise, use as many as are available
		if(nrow(threePlus) < top){
			end3 <- nrow(threePlus)
		} else {
			end3 <- top
		}
		
		top103Plus <- threePlus[1:end3,]
		
		#Perform linear regression on the top 10
		linModel3 <- lm(top103Plus$patternScore ~ seq(1:end3), data = top103Plus)
		
		#Create data frame to hold information
		tempFrame <- data.frame(seq    = inData,
														microhomologyScore = mhScore,
														outOfFrameScore = outOfFrameScore,
														slopeMH3Plus = linModel3$coefficients[2],
														frameShift = fShift,
														stringsAsFactors = FALSE)
	} else {
		tempFrame <- data.frame(seq = "",
														microhomologyScore = 0,
														slopeMH3Plus = 0,
														frameShift = "NA",
														stringsAsFactors = FALSE)
	}
	return(tempFrame)
}

#' calculateBae
#'
#' @param seq The sequence to use for Bae calculation
#' @param cutPosition The index of the nucleotide to the left of the cut site
#' @param weight The length weighting factor to use; default is 20
#'
#' @return
#' @export
#'
#' calculateBae("GTGGCCGACGGGCTCATCACCACGCTCCATTATCCAGCCCCAAAGCGCAA", 25, 20)

#The cutPosition is the index of the base directly to the LEFT of the cut
calculateBae <- function(seq, cutPosition, weight = 20.0){
  seq          <- seq
  left         <- cutPosition
  lengthWeight <- weight
  right        <- nchar(seq) - left + 1
  
  mhDF         <- data.frame(seq     = as.character(),
                             iC      = as.numeric(),
                             ipluskC = as.numeric(),
                             jC      = as.numeric(),
                             jplusk  = as.numeric(),
                             len     = as.numeric(),
                             stringsAsFactors = FALSE)
  
  psDF         <- data.frame(seqDeletion   = as.character(),
                             microhomology = as.character(),
                             delLength     = as.numeric(),
                             patternScore  = as.numeric(),
                             stringsAsFactors = FALSE)
  
  for(k in (left-1):2){
    for(j in right:nchar(seq)){
      for(i in 2:(left-k)){
        if(substring(seq, i, i+k) == substring(seq, j, j+k)){
          length <- j-i
          tempDF <- data.frame(seq     = substring(seq, i, i+k),
                               iC      = i,
                               ipluskC = i+k,
                               jC      = j,
                               jplusk  = j+k,
                               len     = length,
                               stringsAsFactors = FALSE)
          mhDF <- rbind(mhDF, tempDF)
        }
      }
    }
  }
  
  sumScore3    <- 0
  sumScoreNot3 <- 0
  
  if(nrow(mhDF) != 0){
  	for(i in 1:nrow(mhDF)){
  		n          <- 0
  		score3     <- 0
  		scoreNot3  <- 0
  		line       <- mhDF[i,]
  		scrap      <- line$seq
  		leftStart  <- line$iC
  		leftEnd    <- line$ipluskC
  		rightStart <- line$jC
  		rightEnd   <- line$jplusk
  		leng       <- line$len
  		
  		for(j in 1:nrow(mhDF)){
  			if(i != j){
  				lineRef       <- mhDF[j,]
  				leftStartRef  <- lineRef$iC
  				leftEndRef    <- lineRef$ipluskC
  				rightStartRef <- lineRef$jC
  				rightEndRef   <- lineRef$jplusk
  				
  				if((leftStart  >= leftStartRef) &
  					 (leftEnd    <= leftEndRef) &
  					 (rightStart >= rightStartRef) &
  					 (rightEnd   <= rightEndRef)){
  					if(((leftStart - leftStartRef) == (rightStart - rightStartRef)) &
  						 ((leftEnd   - leftEndRef)   == (rightEnd   - rightEndRef))){
  						n <- n + 1
  					}
  				}
  			}
  		}
  		
  		if(n == 0){
  			
  			lengthFactor <- round(1/exp((leng)/lengthWeight), 3)
  			numGC <- str_count(scrap, 'G') + str_count(scrap, 'C')
  			
  			if((leng %% 3) == 0){
  				score3 <- 100 * lengthFactor * ((nchar(scrap) - numGC) + (numGC * 2))
  			} else {
  				scoreNot3 <- 100 * lengthFactor * ((nchar(scrap) - numGC) + (numGC * 2))
  			}
  			
  			tempPsDF <- data.frame(seqDeletion = paste0(substring(seq, 1, leftEnd),
  																									paste(rep('-', leng), collapse = ''),
  																									substring(seq, rightEnd + 1, nchar(seq))),
  														 microhomology = scrap,
  														 delLength     = leng,
  														 patternScore  = (100 * lengthFactor * ((nchar(scrap) - numGC) + (numGC * 2))),
  														 stringsAsFactors = FALSE)
  			
  			psDF <- rbind(psDF, tempPsDF)
  		}
  		
  		sumScore3 <- sumScore3 + score3
  		sumScoreNot3 <- sumScoreNot3 + scoreNot3
  	}
  }
  return(psDF)
}
