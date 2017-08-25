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
calculateBae <- function(seq, cutPosition = 40, weight = 20.0){
	seq          <- seq
	left         <- cutPosition
	lengthWeight <- weight
	right        <- nchar(seq) - left + 1
	
	
	#print(paste0('length of seq = ', nchar(seq)))
	
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
	
	for(k in (left-1):1){
		for(j in right:nchar(seq)){
			for(i in 1:(left-k)){
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
	
	#print(paste0("MH Score: ", sumScore3 + sumScoreNot3))
	#print(paste0("Out-of-frame score: ", (sumScoreNot3*100/(sumScore3 + sumScoreNot3))))
	
	
	return(psDF)
}




