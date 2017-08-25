calculateScores <- function(targets){
	#Define a data frame to hold scoring information
	scoreFrame <- data.frame(id          = as.character(), 
													 seq         = as.character(), 
													 mhScore     = as.numeric(), 
													 oofScore    = as.numeric(), 
													 menthuScore = as.numeric(),
													 stringsAsFactors = FALSE)
	
	#For each target site:
	for(i in targets){
		#Calculate deletion information
		baeFrame <- calculateBae(targets$seq[i])
		
		#Determine Bae microhomology and out-of-frame scores, as well as Menthu score
		tFrame <- data.frame(id          = targets$id[i], 
												 seq         = targets$seq[i],
												 mhScore     = sum(baeFrame$patternScore),
												 oofScore    = sum(baeFrame$patternScore[which(baeFrame$delLength %% 3 != 0)]),
												 menthuScore = calculateMenthu(baeFrame),
												 stringsAsFactors = FALSE)
		

		scoreFrame <- rbind(scoreFrame, tFrame)
	}
	
	return(scoreFrame)
}
