##########Contains functions for dealing with target site identification

#' getPreGenPamList
#'
#' @param pamList 
#'
#' @return
#' @export
#'
#' @examples
#' 

getPreGenPamList <- function(pamList){
	# Pre-defined lists
	cas9Like   <- c('NGG', 'NRG', 'NNNRRT', 'NNGRRT', 'NNGTGA', 'NNAGAAW', 'NNNVRYAC', 'NNNNGMTTT')
	cas12aLike <- c('TTTN', 'TTTV', 'TTN', 'YTN')
	
	# Drop 'redundant' PAMs to save time; retain the broader case
	if("NGG" %in% pamList && "NRG" %in% pamList){
		pamList <- pamList[-which(pamList == "NGG")] 
	}
	
	if("TTTV" %in% pamList && "TTTN" %in% pamList){
		pamList <- pamList[-which(pamList == "TTTN")] 
	}
	
	if("TTN" %in% pamList && "YTN" %in% pamList){
		pamList <- pamList[-which(pamList == "TTN")] 
	}
	
	# Create data frame to hold all this info, and assign cut distances and overhangs
	pamListFrame <- data.frame(pamList = pamList,
														 cutDist = sapply(pamList, function(x) {if(x %in% cas9Like){-3} else if(x %in% cas12aLike){18}}),
														 ohList  = sapply(pamList, function(x) {if(x %in% cas9Like){ 0} else if(x %in% cas12aLike){ 5}}),
														 stringsAsFactors = FALSE)
	
	return(pamListFrame)
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
#' @examples 
#' 
#' @export 
#' 

pamScan <- function(pamList, cutDistList, ohList, targetList, exonList, exonStarts = NULL, findCut = FALSE, type = NULL, wiggle = TRUE, wiggleRoom = 39) {
	require(Biostrings)
	cutDLFrame <- data.frame(Target = pamList, CutDist = cutDistList, stringsAsFactors = FALSE)
	ohLenFrame <- data.frame(Target = pamList, ohLen   = ohList,      stringsAsFactors = FALSE)
	
	if(!is.null(exonStarts)){
		exonMergeFrame <- data.frame(Exon_Num = exonList, exonStarts = exonStarts, stringsAsFactors = FALSE)

	}
	
	#Convert IUPAC PAMs to regular expressions
	pamPosRegs <- sapply(strsplit(pamList, ""),                         function(x) paste0("[", paste(Biostrings::IUPAC_CODE_MAP[x], collapse = "]["), "]"))
	pamNegRegs <- sapply(strsplit(reverseComplement.list(pamList), ""), function(x) paste0("[", paste(Biostrings::IUPAC_CODE_MAP[x], collapse = "]["), "]"))
	
	#Add lookahead to PAMs
	pamPosRegs <- paste0("(?=", pamPosRegs, ")")
	pamNegRegs <- paste0("(?=", pamNegRegs, ")")
	
	#Find all PAM matches in all sequences
	matchesPos <- sapply(pamPosRegs, function(x) sapply(targetList, function(y) gregexpr(x, y, perl = TRUE)))
	matchesNeg <- sapply(pamNegRegs, function(x) sapply(targetList, function(y) gregexpr(x, y, perl = TRUE)))
	
	#If there is only a single exon and a single PAM
	if(class(matchesPos) == "list" && length(pamList) == 1){
		#Extract the loci of the matches
		matchesPos <- sapply(1:length(matchesPos), function(x) matchesPos[[x]][1:length(matchesPos[[x]])])
		matchesNeg <- sapply(1:length(matchesNeg), function(x) matchesNeg[[x]][1:length(matchesNeg[[x]])])
		
		#Format the results matches
		pamFramePos <- data.frame(Target           = rep(pamList,        length(matchesPos)),
															Orientation      = rep("forward",      length(matchesPos)),
															Exon_Num         = rep(exonList,       length(matchesPos)),
															Sites            = matchesPos,
															MatchLength      = rep(nchar(pamList), length(matchesPos)),
															CutIndex         = rep(0,              length(matchesPos)),
															ohLength         = rep(1,              length(matchesPos)),
															contextCondition = rep(FALSE,          length(matchesPos)),
															stringsAsFactors = FALSE)
		
		#Format the results matches
		pamFrameNeg <- data.frame(Target           = rep(pamList,        length(matchesNeg)),
															Orientation      = rep("complement",   length(matchesNeg)),
															Exon_Num         = rep(exonList,       length(matchesNeg)),
															Sites            = matchesNeg + 1,
															MatchLength      = rep(nchar(pamList), length(matchesNeg)),
															CutIndex         = rep(0,              length(matchesNeg)),
															ohLength         = rep(2,              length(matchesNeg)),
															contextCondition = rep(FALSE,          length(matchesNeg)),
															stringsAsFactors = FALSE)
		
	} else if(class(matchesPos) == "list" && length(pamList) > 1){	
		#Extract the loci of the matches
		matchesPos <- sapply(1:length(matchesPos), function(x) matchesPos[[x]][1:length(matchesPos[[x]])])
		matchesNeg <- sapply(1:length(matchesNeg), function(x) matchesNeg[[x]][1:length(matchesNeg[[x]])])
		
		#Rename the lists
		names(matchesPos) <- pamList
		names(matchesNeg) <- pamList
		
		#Stack the lists
		posStack <- stack(matchesPos)
		negStack <- stack(matchesNeg)
		
		#Format the results matches
		pamFramePos <- data.frame(Target           =        as.character(posStack[ , 2]),
															Orientation      = rep("forward", nrow(posStack)),
															Exon_Num         = rep(exonList,  nrow(posStack)),
															Sites            = posStack[ , 1],
															MatchLength      = nchar(as.character(posStack[ , 2])),
															CutIndex         = rep(0,         nrow(posStack)),
															ohLength         = rep(3,         nrow(posStack)),
															contextCondition = rep(FALSE,     nrow(posStack)),
															stringsAsFactors = FALSE)
		
		pamFrameNeg <- data.frame(Target           =           as.character(negStack[ , 2]),
															Orientation      = rep("complement", nrow(negStack)),
															Exon_Num         = rep(exonList,     nrow(negStack)),
															Sites            = negStack[ , 1] + 1,
															MatchLength      = nchar(as.character(negStack[ , 2])),
															CutIndex         = rep(0,            nrow(negStack)),
															ohLength         = rep(4,            nrow(negStack)),
															contextCondition = rep(FALSE,        nrow(negStack)),
															stringsAsFactors = FALSE)
		
	} else {
		
	#Rename resulting matrix columns (PAMs)
	colnames(matchesPos) <- pamList	
	colnames(matchesNeg) <- pamList
	
	#Rename resulting matrix rows (target sequences)
	row.names(matchesPos) <- exonList
	row.names(matchesNeg) <- exonList
	
	posStack <- stack(matchesPos)
	negStack <- stack(matchesNeg)
	
	#Create data frame to hold all PAM info
	pamFramePos <- data.frame(Target           = as.character(posStack[ , 2]),
														Orientation      = rep("forward", nrow(posStack)),
														Exon_Num         = posStack[ , 1],
														Sites            = posStack[ , 4],
														MatchLength      = nchar(posStack[ , 2]),
														CutIndex         = rep(0,     nrow(posStack)),
														ohLength         = rep(4,     nrow(posStack)),
														contextCondition = rep(FALSE, nrow(posStack)),
														stringsAsFactors = FALSE)
	
	pamFrameNeg <- data.frame(Target           = as.character(negStack[ , 2]),
														Orientation      = rep("complement", nrow(negStack)),
														Exon_Num         = negStack[ , 1],
														Sites            = negStack[ , 4] + 1,
														MatchLength      = nchar(negStack[ , 2]),
														CutIndex         = rep(0,     nrow(negStack)),
														ohLength         = rep(4,     nrow(negStack)),
														contextCondition = rep(FALSE, nrow(negStack)),
														stringsAsFactors = FALSE)
	}

	# Drop -1 returns for not found instances
	pamFramePos <- pamFramePos[which(pamFramePos$Sites > -1), ]
	pamFrameNeg <- pamFrameNeg[which(pamFrameNeg$Sites > -1), ]
	
	# Add cut distance to each row based on PAM
	pamFramePos <- plyr::join(pamFramePos, cutDLFrame, by = "Target")
	pamFrameNeg <- plyr::join(pamFrameNeg, cutDLFrame, by = "Target")
	
	# Add overhang length to each row based on PAM
	pamFramePos <- plyr::join(pamFramePos, ohLenFrame, by = "Target")
	pamFrameNeg <- plyr::join(pamFrameNeg, ohLenFrame, by = "Target")
	
	if(!is.null(exonStarts)){
		# Add exon starts to each exon
		pamFramePos <- plyr::join(pamFramePos, exonMergeFrame, by = "Exon_Num")
		pamFrameNeg <- plyr::join(pamFrameNeg, exonMergeFrame, by = "Exon_Num")
	}
	
	# If we are not interested in finding the cut site
	if(!findCut){
		pamFramePos$CutIndex <- pamFramePos$Sites
		pamFrameNeg$CutIndex <- pamFrameNeg$Sites
	
	# Find the cut site from the match index and other information
	} else {
		
		posFuncAdj <- function(x) {
			  pamFramePos$Sites[x]   + 
				pamFramePos$CutDist[x] +
				if(pamFramePos$CutDist[x] < 0){
					-1
				} else {
					pamFramePos$MatchLength[x] - 1
				}
		}
		
		negFuncAdj <- function(x){
		  	pamFrameNeg$Sites[x]   - 
				pamFrameNeg$CutDist[x] +
				
				if(pamFrameNeg$CutDist[x] < 0){
					pamFrameNeg$MatchLength[x] - 2
				} else {
					-2
				}
		}
		
		pamFramePos$CutIndex <- sapply(1:nrow(pamFramePos), posFuncAdj)
		pamFrameNeg$CutIndex <- sapply(1:nrow(pamFrameNeg), negFuncAdj)	                      			 	

	}
	
	# Merge the data frames
	pamFrame <- rbind(pamFramePos, pamFrameNeg)
	

	# Correct positions by localization in gene and not exon
	if(!is.null(exonStarts)){
		if(wiggle){
			
			pamFrame$CutIndex <- sapply(1:nrow(pamFrame), 
																	function(x) pamFrame$CutIndex[x] + pamFrame$exonStarts[x] - (if(pamFrame$exonStarts[x] - wiggleRoom > 0){wiggleRoom} else {0}) - 1)
		} 
	}
	
	return(pamFrame)
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

talPal <- function(targetList, findCut = TRUE, wiggle = TRUE, wiggleRoom = 39, range = TRUE, 
									 armin = 15, armax = 18, spamin = 14, spamax = 16, exonList, exonStarts = NULL, Exon_Num = NULL) {
	require(stringr)
	require(plyr)
	
	targetList <- as.character(targetList) # Transforms DNAString to char
	
	if(!is.null(exonStarts)){
		exonMergeFrame <- data.frame(Exon_Num = exonList, exonStarts = exonStarts, stringsAsFactors = FALSE)
		
	}
	
	# Variable definition and initialization
	numSites <- length(targetList) # Number of sites to target
	talFrame <- data.frame(talSites         = numeric(numSites), 
												 talSeqs          = character(numSites), 
												 stringsAsFactors = FALSE)
	
	# Calculates range of allowable spacer and arm lengths
	if (range == TRUE) {
		#Generate all possible search patterns given the TALEN parameters
		patternLists <- getTalPatterns(armin, armax, spamin, spamax) 
		
		patternPos   <- patternLists[[1]] #Search patterns for forward strand
		patternNeg   <- patternLists[[2]] #Search patterns for complementary strand
		armSpaList   <- patternLists[[3]] #Arm/spacer/arm length labels
		
		#Find all TALEN matches in all sequences
		matchesPosT <- sapply(patternPos, function(x) sapply(targetList, function(y) gregexpr(x, y, perl = TRUE)))
		#matchesNegT <- sapply(patternNeg, function(x) sapply(targetList, function(y) gregexpr(x, y, perl = TRUE)))
		
		if(class(matchesPosT) == "list"){
			#Extract the matches
			matchesPosT <- sapply(1:length(matchesPosT), function(x) matchesPosT[[x]][1:length(matchesPosT[[x]])])
			#matchesNegT <- sapply(1:length(matchesNegT), function(x) matchesNegT[[x]][1:length(matchesNegT[[x]])])
			
			#Rename the lists
			names(matchesPosT) <- armSpaList
			#names(matchesNegT) <- armSpaList
			
			#Stack the lists
			posStackT <- stack(matchesPosT)
			#negStackT <- stack(matchesNegT)
			
			#Get the length of the arms+spacer for the TALEN
			matchLengthsPos <- sapply(1:nrow(posStackT), function(x) sum(as.numeric(unlist(strsplit(as.character(posStackT[x, 2]), "/")))))
			#matchLengthsNeg <- sapply(1:nrow(negStackT), function(x) sum(as.numeric(unlist(strsplit(as.character(negStackT[x, 2]), "/")))))
			
			#Format the results matches
		  talFramePos <- data.frame(Target           = as.character(posStackT[ , 2]),
																Orientation      = rep("forward",  nrow(posStackT)),
																Exon_Num         = rep(1,          nrow(posStackT)),
																Sites            = posStackT[ , 1],
																MatchLength      = nchar(as.character(posStackT[ , 2])),
																CutIndex         = rep(0,        nrow(posStackT)),
																contextCondition = rep(FALSE,    nrow(posStackT)),
																stringsAsFactors = FALSE)
			
			# talFrameNeg <- data.frame(Target           = as.character(negStackT[ , 2]),
			# 													Orientation      = rep("complement", nrow(negStackT)),
			# 													Exon_Num         = rep(1,            nrow(negStackT)),
			# 													Sites            = negStackT[ , 1],
			# 													MatchLength      = nchar(as.character(negStackT[ , 2])),
			# 													CutIndex         = rep(0,        nrow(negStackT)),
			# 													contextCondition = rep(FALSE,    nrow(negStackT)),
			# 													stringsAsFactors = FALSE)
			
		} else {
			#Rename resulting matrix columns (PAMs)
			colnames(matchesPosT) <- armSpaList	
			#colnames(matchesNegT) <- armSpaList
			
			#Rename resulting matrix rows (target sequences)
			row.names(matchesPosT) <- exonList
			#row.names(matchesNegT) <- exonList
			
			#Unstack matrix into nice data frame
			posStackT <- stack(matchesPosT)
			#negStackT <- stack(matchesNegT)
			
			#Get the length of the arms+spacer for the TALEN
			matchLengthsPos <- sapply(1:nrow(posStackT), function(x) sum(as.numeric(unlist(strsplit(as.character(posStackT[x, 2]), "/")))))
			#matchLengthsNeg <- sapply(1:nrow(negStackT), function(x) sum(as.numeric(unlist(strsplit(as.character(negStackT[x, 2]), "/")))))
			
			#Create data frame to hold all TALEN info
			talFramePos <- data.frame(Target           = as.character(posStackT[ , 2]),
																Orientation      = rep("forward", nrow(posStackT)),
																Exon_Num         = posStackT[ , 1],
																Sites            = posStackT[ , 4],
																MatchLength      = matchLengthsPos,
																CutIndex         = rep(0,        nrow(posStackT)),
																contextCondition = rep(FALSE,    nrow(posStackT)),
																stringsAsFactors = FALSE)
			
			# talFrameNeg <- data.frame(Target           = as.character(negStackT[ , 2]),
			# 													Orientation      = rep("complement", nrow(negStackT)),
			# 													Exon_Num         = negStackT[ , 1],
			# 													Sites            = negStackT[ , 4] + 1,
			# 													MatchLength      = matchLengthsNeg,
			# 													CutIndex         = rep(0,            nrow(negStackT)),
			# 													contextCondition = rep(FALSE,        nrow(negStackT)),
			# 													stringsAsFactors = FALSE)
			
		}

		#Drop -1 returns for not found instances
		talFramePos <- talFramePos[which(talFramePos$Sites > -1), ]
		#talFrameNeg <- talFrameNeg[which(talFrameNeg$Sites > -1), ]
		
		if(!findCut){
			talFramePos$CutIndex <- talFramePos$Sites
			#talFrameNeg$CutIndex <- talFrameNeg$Sites
			
		} else {
			storePos <- sapply(1:length(talFramePos$Target), function (x) list(as.numeric(unlist(strsplit(as.character(talFramePos$Target[x]), "/")))))
			#storeNeg <- sapply(1:length(talFrameNeg$Target), function (x) list(as.numeric(unlist(strsplit(as.character(talFrameNeg$Target[x]), "/")))))
			
			#Get left spacer arm + half the spacer
			cutModPos <- talFramePos$Sites + sapply(1:length(storePos), function(x) as.numeric(storePos[[x]][1]) + (as.numeric(storePos[[x]][2]) / 2))
			#cutModNeg <- talFrameNeg$Sites + sapply(1:length(storeNeg), function(x) as.numeric(storeNeg[[x]][1]) + (as.numeric(storeNeg[[x]][2]) / 2))
			
			#Add cut distance to each row based on PAM
			talFramePos$CutIndex <- cutModPos
			#talFrameNeg$CutIndex <- cutModNeg
		}
		
		if(!is.null(exonStarts)){
			#Add exon starts to each exon
			talFramePos <- plyr::join(talFramePos, exonMergeFrame, by = "Exon_Num")
			#talFrameNeg <- plyr::join(talFrameNeg, exonMergeFrame, by = "Exon_Num")
		}
		
		#Concatenate the data frames
		#talFrame <- rbind(talFramePos, talFrameNeg)
		talFrame  <- talFramePos
		
		# Correct positions by localization in gene and not exon, taking into account "wiggle room"
		if (!is.null(exonStarts)) {
			if(wiggle){
				talFrame$CutIndex <- sapply(1:nrow(talFrame), 
																		function(x) talFrame$CutIndex[x] + talFrame$exonStarts[x] - (if(talFrame$exonStarts[x] - wiggleRoom > 0){wiggleRoom} else {0}) - 1)
			} else {
				talFrame$CutIndex <- sapply(1:nrow(talFrame), 
																		function(x) talFrame$CutIndex[x] + talFrame$exonStarts[x])
			}
		}
	}
	
	return(talFrame)
}

getTalPatterns <- function(armin, armax, spamin, spamax){
	patListPos <- NULL
	patListNeg <- NULL
	patList    <- NULL
	
	patListPos <- sapply(spamin:spamax, function(x) sapply(armin:armax, function(y) sapply(armin:armax, function(z) paste0("(?=(T[ACGT]{", x + y + z - 2, "}A))"))))
	patListNeg <- sapply(spamin:spamax, function(x) sapply(armin:armax, function(y) sapply(armin:armax, function(z) paste0("(?=(A[ACGT]{", x + y + z - 2, "}T))"))))
	patList    <- sapply(spamin:spamax, function(x) sapply(armin:armax, function(y) sapply(armin:armax, function(z) paste0(y, "/", x, "/", z))))
	
	if(class(patListPos) != "character"){
		if(ncol(patListPos > 2)){
			patListPos <- patListPos[, -2]
			patListNeg <- patListNeg[, -2]
			patList    <-    patList[, -2]
		}
	}

	patListPos <- unlist(as.list(patListPos))
	patListNeg <- unlist(as.list(patListNeg))
	patList    <- unlist(as.list(patList))
	
	return(list(patListPos, patListNeg, patList))
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
		tempL  <- NULL
		tempR  <- NULL
		spacer <- NULL
		talL   <- NULL
		talR   <- NULL
		spa    <- NULL
		specs  <- NULL
		
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
