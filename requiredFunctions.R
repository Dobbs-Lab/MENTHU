####Required Function#######

calculateMENTHUGeneSeq <- function(casList, geneSeq, threshold, exons, progress){
	pamList <- casList
	#pamList <- "NGG"
	#print(pamList)
	#exon <- getExon("GU179289.1", percent = 10)
	#exonInfo <- exon[[1]]
	#exonSeq <- exon[[2]]
	#genSeq <- exon[[3]]
	
	#	pamList <- list("NGG","NGCG", "NGAG", "NGAN", "NGNG", "NNGRRT", "NNGRRN", "NNNNGATT", "NNAGAAW", "NAAAAC")
	
	#pamSites <- pamScan(pamList,DNAStringSet(input$geneSeq),findCut = TRUE,type = "cas9")
	
	#Subset geneSeq input to only search exons for PAMs
	#exonSeqs <- geneSeq
	progress$inc(0.1, detail = "Scanning for target sites...")
	pamSites <- pamScan(pamList,DNAStringSet(geneSeq),findCut = TRUE,type = "cas9")	
	
	#print(pamSites)
	#test <- pamSites[[1]][[1]][[1]][[1]]
	
	# print(test)
	#print(window(genSeq,test,80))
	
	menthuFrame <- data.frame(targetSequence = as.character(), menthuScore = as.numeric(), toolType = as.character(), strand = as.character(), exon = as.character(), location = as.integer())
	#PAM LEVEL
	for(i in 1:length(pamSites[])){
		progress$inc(1/length(pamSites[]), detail = paste("Processing ", names(pamSites[i])))
		toolTypeI <- names(pamSites)[i]
		
		#STRAND LEVEL
		#for(j in 1:length(pamSites[[i]][])){
		
		#EXON LEVEL
		for(k in 1:length(pamSites[[i]][[1]][])){
			exonNum <- names(pamSites[[i]][[1]])[k]
			#SITE LEVEL
			for(l in 1:length(pamSites[[i]][[1]][[k]][])){
				print(l)
				#If there's enough context...
				if(!is.na(pamSites[[i]][[1]][[k]][[l]]) && !is.null(pamSites[[i]][[1]][[k]][[l]])){
					if((pamSites[[i]][[1]][[k]][[l]] > 43) && ((pamSites[[i]][[1]][[k]][[l]] + 40 + nchar(toolTypeI)) < nchar(geneSeq))){
						#context <- window(DNAString(input$geneSeq), pamSites[[i]][[1]][[k]][[l]], 80)
						context <- window(DNAString(geneSeq), pamSites[[i]][[1]][[k]][[l]], 80)
						slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
						print(abs(slopeFrame$slopeMH3Plus))
						
						#if(abs(slopeFrame$slopeMH3Plus) > input$threshold){
						if(abs(slopeFrame$slopeMH3Plus) >= threshold){
							crispr <- substr(slopeFrame$seq, 24, (43 + nchar(toolTypeI)))
							
							formFrame  <- data.frame(targetSequence = crispr, menthuScore = round(abs(slopeFrame$slopeMH3Plus), digits = 2), toolType = toolTypeI, strand = "forward", exon = exonNum, location = as.integer(pamSites[[i]][[1]][[k]][[l]], digits = 0))
							menthuFrame <- rbind(menthuFrame, formFrame)
						}
					}
				}
			}
		}
		#}
	}
	print(menthuFrame)
	#colnames(menthuFrame) <- c("Target Sequence", "Score", "Tool", "Strand", "Exon", "Location")
	return(menthuFrame)
}


####Exon Handler for Custom Exon Input####
exonHandler <- function(exonRHandsonTable){
	exonTable <- hot_to_r(exonRHandsonTable)
	exons <- exonTable[apply(exonTable, MARGIN = 1, function(x) any(x > 0)),]
	
	return(exons)
}

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
pamScan <- function(pamList, targetList, exonStarts = NULL, findCut = FALSE, type = NULL) {
	
	require(Biostrings)
	
	# Variable definition and initialization
	numPAM <- length(pamList) # Number of PAMs in pamList
	pamSites <- vector("list", numPAM)
	numTarget <- length(targetList) # Number of DNA targets to screen
	
	for (i in 1:length(pamList)) { # Iterate for each PAM
		
		# Variable initialization
		ppam <- DNAString(pamList[[i]]) # Transform char to DNAString
		npam <- reverseComplement(ppam) # Find PAM sequence in - DNA strand
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
				switch(type,
							 cas9 = {pShift = -4}
				)
				switch(type,
							 cas9 <- {nShift = 5}
				)
				for (k in 1:length(pdf$start)) {
					ppamSites[k] <- pdf$start[k] + pShift
				}
				for (k in 1:length(ndf$start)) {
					npamSites[k] <- ndf$start[k] + nShift
				}
			}
			
			
			# Add ppamStes and npamSites to memory
			pPAM[[j]] <- ppamSites
			nPAM[[j]] <- npamSites
			
			# Correct positions by localization in gene and not exon
			if (is.null(exonStarts) == FALSE) {
				pPAM[[j]] <- pPAM[[j]] + exonStarts[j]
				nPAM[[j]] <- nPAM[[j]] + exonStarts[j]
			}
		}
		
		# Name exons
		names(pPAM)[1:numTarget] <- paste0("Exon ", 1:(numTarget))
		names(nPAM)[1:numTarget] <- paste0("Exon ", 1:(numTarget))
		
		# Name PAMs by strand
		PAMs <-list(pPAM, nPAM)
		names(PAMs)[1:2] <- paste0(pamList[[i]], c(": Positive strand", ": Negative strand"))
		
		# Name exons by PAM
		pamSites[[i]] <- PAMs
		names(pamSites)[i] <- pamList[[i]]
	}
	return(pamSites)
}


#' getExon
#' 
#' Function that fetches exon data from a GenBank accesion ID
#' 
#' @param accession GenBank accession ID
#' @param percent Percentage of exons from 5' to 3' to fetch
#' 
#' @result List of exon start and end site and DNA sequences
#' 
#' @examples 
#' 
#' @export 
#' 
window <- function(sequence, position, winSize = 80) {
	
	if(position < winSize/2){
		return("")
	} else {
		return(sequence[(position-((winSize + winSize%%2)/2 - 1)):(position+((winSize + winSize%%2)/2))])
	}
}

#' calculateSlopeCompetition
#'
#' @param inData A table of two columns - first is geneID, second is wildtype sequence
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight A weight factor for caculating pattern score; default is 20, as in Bae paper
#' @param top The number of entries to consider when generating the linear regression
#'
#' @return


#For each sequence in the list, calculate the slope competition
calculateSlopeCompetition <- function(inData, cutSite = -1, weight = 20, top = 10){
	#Create data frame to hold information about targets
	#targetDF <- data.frame(seq                = as.character(),
	#											 microhomologyScore = as.numeric(),
	#											 outOfFrameScore    = as.numeric(),
												 #slopeMH2Plus       = as.numeric(),
	#											 slopeMH3Plus       = as.numeric(),
	#											 stringsAsFactors = FALSE)
	
	lengthWeight <- weight
	
	#for(a in 1:nrow(inData)){
		if(cutSite == -1){
			#If the cut site isn't specified, assume it is in the middle of the input
			cutSite <- floor(nchar(inData)/2)
		}
		
		#Get the pattern scores for every deletion pattern in the wildtype sequence
		patternScoreDF <- calculateBae(inData, cutSite, lengthWeight)
		patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore),]
		
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
		
		threePlus <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3),]
		if(nrow(threePlus) < top){
			end3 <- nrow(threePlus)
		} else {
			end3 <- top
		}
		top103Plus <- threePlus[1:end3,]
		
		#if(nrow(patternScoreDF) < top){
		#  end2 <- nrow(patternScoreDF)
		#} else {
		#  end2 <- top
		#}
		
		#Perform linear regression on the top 10
		linModel3 <- lm(top103Plus$patternScore ~ seq(1:end3), data = top103Plus)
		#linModel2 <- lm(patternScoreDF$patternScore[1:end2] ~ seq(1:end2), data = patternScoreDF)
		#Create data frame to hold information
		tempFrame <- data.frame(seq    = inData,
														microhomologyScore = mhScore,
														outOfFrameScore = outOfFrameScore,
														#slopeMH2Plus = linModel2$coefficients[2],
														slopeMH3Plus = linModel3$coefficients[2],
														stringsAsFactors = FALSE)
		
		#targetDF <- rbind(targetDF, tempFrame)
	#}
	
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




