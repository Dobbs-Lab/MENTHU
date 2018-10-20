####Required Functions#######

#' calculateMENTHUGeneSeq
#'
#' @param casList 
#' @param wiggle 
#' @param wiggleRoom 
#' @param geneSeq 
#' @param threshold 
#' @param exonDF 
#' @param progress 
#' @param armin 
#' @param armax 
#' @param spamin 
#' @param spamax 
#'
#' @return
#' @export
#'
#' @examples

#calculateMENTHUGeneSeq <- function(casList, cutDistList, wiggle = TRUE, wiggleRoom = 39, geneSeq, threshold, talFlag,
#																	 exonDF, progress, armin, armax, spamin, spamax, version){
#calculateMENTHUGeneSeq <- function(casList, cutDistList, wiggle = TRUE, wiggleRoom = 39, geneSeq, threshold, talFlag,
#																	 exonDF, progress, armin, armax, spamin, spamax){
calculateMENTHUGeneSeq <- function(casList, cutDistList, wiggle = TRUE, wiggleRoom = 39, geneSeq, exonDF, progress){
	require(Biostrings)
	require(plyr)
	
	# Set variables
	talFlag   <- FALSE
	version   <- 2
	noPamFlag <- FALSE
	
	geneSeq <- toupper(geneSeq)
	
	#Rename the input casList to pamList, because I'm currently too lazy to change variable names
	pamList <- casList
	
	#If both NGG and NRG are selected in the PAM list, just search for NRG to save time
	if("NGG" %in% pamList & "NRG" %in% pamList){
		pamList     <- pamList[-1]
		cutDistList <- cutDistList[-1]
	}
	
	########Subset geneSeq to exon sequences############
	# Deal with the case in which exons are specified and extra context should be included to examine sites where the gRNA would run off the exon
	if((wiggle == TRUE) && (class(exonDF) == "data.frame")){
		# Ensure extra context doesn't run off the end of the whole sequence
		exStart <- sapply(1:nrow(exonDF), function(x) if((exonDF$exonStart[x] - wiggleRoom) < 1){1}                           else {exonDF$exonStart[x] - wiggleRoom})
		
		# Ensure extra context doesn't run off the end of the whole sequence
		exEnd   <- sapply(1:nrow(exonDF), function(x) if((exonDF$exonEnd[x]   + wiggleRoom) > nchar(geneSeq)){nchar(geneSeq)} else {exonDF$exonEnd[x]   + wiggleRoom})
		
		# Create new exonSeqs with 'fixed' context
		exonSeqs <- substring(rep(toupper(geneSeq), length(exStart)), exStart, exEnd)
		
		# Deal with the case in which exons are specified but extra context should not be included
	} else if((!wiggle) && (class(exonDF) == "data.frame")){
		exonSeqs <- substring(rep(toupper(as.character(geneSeq)), nrow(exonDF)), 
													exonDF$exonStart, 
													exonDF$exonEnd)
		
	} else {
		# In the case where exons are not specified, treat the whole gene sequence as a single exon
		exonSeqs <- geneSeq
		
	}
	
	# If the user is using Cas:
	if(length(pamList) > 0){
		
		# Update progress bar
		progress$inc(0.01, detail = "Scanning for CRISPR target sites...")
		
		# If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
		if(class(exonDF) == "data.frame"){
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeqs, 
													exonList   = exonDF$Exon_Num, 
													exonStarts = exonDF$exonStart,
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wiggleRoom = wiggleRoom)
			
		} else {
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeqs,
													exonList   = 1,
													exonStarts = 1, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wiggleRoom = wiggleRoom)
		}
		
		# Count the number of Cas target sites
		siteCount <- nrow(pamSites)
		
		# Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
		# Update progress bar
		progress$inc(0.01, detail = "Pre-processing CRISPR target sites...")
		
		if(class(exonDF) == "data.frame"){
			pamSites <- unique(suppressMessages(plyr::join(pamSites, exonDF, by = 'Exon_Num')))
			
		} else {
			tempExonDF <- data.frame(Exon_Num  = 1,
															 exonStart = 1,
															 exonEnd   = nchar(exonSeqs),
															 stringsAsFactors = FALSE)
			
			pamSites <- unique(suppressMessages(plyr::join(pamSites, tempExonDF, by = 'Exon_Num')))
			
		}
		
		# Drop target sites where the cut site is not within the exon boundaries
		keep     <- sapply(1:nrow(pamSites), function(x) pamSites$CutIndex[x] %in% seq(from = pamSites$exonStart[x], to = pamSites$exonEnd[x], by = 1))
		pamSites <- pamSites[keep, ]
		
		# Identify sites with enough sequence context to do calculations
		pamSites$contextCondition[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		pamSites      <- pamSites[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		# Count the number of targets satisfying context condition
		siteCountC <- nrow(pamSites)
		
		if(nrow(pamSites) >= 1){
			# Get the sequence context surrounding the cut site
			context <- unlist(lapply(1:nrow(pamSites), function (x) substr(geneSeq, pamSites$CutIndex[x] - 39, pamSites$CutIndex[x] + 40)))
			
			# Set the sequence context in the frame
			pamSites$seq = context
			
		} else {
			noPamFlag <- TRUE
			
		}
	} else {
		#If the user is NOT using CRISPR/Cas system, set pamFlag to FALSE
		pamSites <- 0
		pamFlag  <- FALSE
		
	}
	
	if(!noPamFlag){
		# Set a flag to be true if there are TALEN inputs
		#talFlag <- armin != "" && armax != "" && spamin != "" && spamax != ""
		
		# If there are TALEN inputs
		if(talFlag){
			
			# Set the range flag to true
			rFlag <- TRUE
			
			# If there are exon inputs
			if(class(exonDF) == "data.frame"){
				# Set all exon starts to the exon starts in the input frame
				# Submit talen info to talPal
				talSites <- talPal(exonSeqs,
													 findCut    = TRUE,
													 wiggle     = TRUE,
													 wiggleRoom    = 39,
													 range      = rFlag, 
													 armin      = armin, 
													 armax      = armax, 
													 spamin     = spamin, 
													 spamax     = spamax, 
													 exonStarts = exonDF$exonStart,
													 exonList   = exonDF$Exon_Num)
				
			} else {
				# If there are no exon inputs, make exon start null
				# Submit talen info to talPal
				talSites <- talPal(exonSeqs,
													 findCut    = TRUE,
													 wiggle     = TRUE,
													 wiggleRoom    = 39,
													 range      = rFlag, 
													 armin      = armin, 
													 armax      = armax, 
													 spamin     = spamin, 
													 spamax     = spamax, 
													 exonStarts = NULL,
													 exonList   = 1)
			}
			
			# Update progress bar
			progress$inc(0.01, detail = "Pre-processing TALEN target sites...")
			
			if(class(exonDF) == "data.frame"){
				talSites <- unique(suppressMessages(plyr::join(talSites, exonDF, by = 'Exon_Num')))
				
			} else {
				tempExonDF <- data.frame(Exon_Num  = 1,
																 exonStart = 1,
																 exonEnd   = nchar(exonSeqs),
																 stringsAsFactors = FALSE)
				
				talSites <- unique(suppressMessages(plyr::join(talSites, tempExonDF, by = 'Exon_Num')))
				
			}
			
			# Drop target sites where the cut site is not within the exon boundaries
			keepT     <- sapply(1:nrow(talSites), function(x) talSites$CutIndex[x] %in% seq(from = talSites$exonStart[x], to = talSites$exonEnd[x], by = 1))
			talSites  <- talSites[keepT, ]
			
			# Identify sites with enough sequence context to do calculations
			talSites$contextCondition[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
			talSites      <- talSites[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq))), ]
			
			# Get the sequence context surrounding the cut site
			contextT <- unlist(lapply(1:nrow(talSites), function (x) substr(geneSeq, talSites$CutIndex[x] - 39, talSites$CutIndex[x] + 40)))
			
			# Set the sequence context in the frame
			talSites$seq = contextT
			
		} else {
			# If TALENs are not used, set talSites to 0
			talSites <- 0
			
		}
		
		# If there is no exon data input, create a dummy data frame to have the 'exon' be the entire gene sequence, starting on nt 1
		if(class(exonDF) != "data.frame"){
			exonDF <- data.frame(exonStart        = 1, 
													 exonEnd          = nchar(geneSeq), 
													 stringsAsFactors = FALSE)
		}
		
		if(version == 1){
			if(pamFlag){
				# Update progress bar
				progress$inc(0.01, detail = "Calculating MENTHUv1.0 scores for CRISPR sites...")
				
				slopeFrameFunc <- function(x){
					# Increment progress
					if(talFlag){
						progress$inc(1 / (nrow(pamSites) + nrow(talSites)))
						
					} else {
						progress$inc(1 / nrow(pamSites))
						
					}
					
					#Return calculation
					return(calculateSlopeCompetition(as.character(x), cutSite = 40, weight = 20, top = 10))
				}
				
				# Calculate slope competition on all the context
				slopeFrame <- as.data.frame(matrix(unlist(sapply(context, slopeFrameFunc)), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
				
				# Update progress bar
				progress$inc(0.01, detail = "Formatting CRISPR site results...")
				
				# Clean up the resulting data frame
				colnames(slopeFrame) <- c("seq", "microhomology_score", "OOF_Score", "slopeMH3Plus", "frameShift", "topDel")
				rownames(slopeFrame) <- c()
				
				# Merge slope frame and pamSites
				pamSites <- unique(suppressMessages(plyr::join(pamSites, slopeFrame)))
				
				# Generate the crispr target sequences
				pamSites$crispr <- sapply(1:nrow(pamSites), 
																	function(x) substring((if(pamSites$Orientation[x] == "forward"){pamSites$seq[x]} else {reverseComplement(pamSites$seq[x])}),
																												40 - pamSites$CutDist[x] - 19,
																												40 - pamSites$CutDist[x] + nchar(pamSites$Target[x])))
				# Clean the data frame
				row.names(pamSites) <- c()
				
				# Create a new data frame of the results
				pamFormFrame <- data.frame(Target_Sequence = pamSites$crispr,
																	 MENTHU_Score    = round(abs(as.numeric(pamSites$slopeMH3Plus)), digits = 2),
																	 Frame_Shift     = pamSites$frameShift,
																	 Tool_Type       = pamSites$Target,
																	 Strand          = pamSites$Orientation,
																	 Exon_ID         = pamSites$Exon_Num,
																	 Cut_Location    = pamSites$CutIndex,
																	 Top_Deletion    = pamSites$topDel,
																	 stringsAsFactors = FALSE)
			}
			
			if(talFlag){
				# Update progress bar
				progress$inc(0.01, detail = "Calculating MENTHUv1.0 scores for TALEN sites...")
				
				slopeFrameTFunc <- function(x){
					# Increment progress
					if(pamFlag){
						progress$inc(1 / (nrow(talSites) + nrow(pamSites)))
						
					} else {
						progress$inc(1 / nrow(talSites))
						
					}
					
					# Return calculation
					return(calculateSlopeCompetition(as.character(x), weight = 20, top = 10))
				}
				
				# Calculate slope competition on all the context
				slopeFrameT <- as.data.frame(matrix(unlist(sapply(contextT, slopeFrameTFunc)), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
				
				# Update progress bar
				progress$inc(0.01, detail = "Formatting TALEN site results...")
				
				# Clean up the resulting data frame
				colnames(slopeFrameT) <- c("seq", "microhomology_score", "OOF_Score", "slopeMH3Plus", "frameShift", "topDel")
				rownames(slopeFrameT) <- c()
				
				# Merge slope frame and pamSites
				talSites <- unique(suppressMessages(plyr::join(talSites, slopeFrameT)))
				
				talenGenFunc <- function(talRow){
					dim  <- unlist(strsplit(talRow$Target, "/"))
					arm1 <- as.numeric(dim[1])
					spa1 <- as.numeric(dim[2])
					arm2 <- as.numeric(dim[3])
					
					armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
					spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
					armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
					
					return(paste0("<strong>", armL, "</strong>", spac, "<strong>", armR, "</strong>"))
				}
				
				# Generate the crispr target sequences
				talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
				
				# Clean the data frame
				row.names(talSites) <- c()
				
				# Create a new data frame of the results
				talFormFrame <- data.frame(Target_Sequence = talSites$talen,
																	 MENTHU_Score    = round(abs(as.numeric(talSites$slopeMH3Plus)), digits = 2),
																	 Frame_Shift     = talSites$frameShift,
																	 Tool_Type       = talSites$Target,
																	 Strand          = talSites$Orientation,
																	 Exon_ID         = talSites$Exon_Num,
																	 Cut_Location    = talSites$CutIndex,
																	 Top_Deletion    = talSites$topDel,
																	 stringsAsFactors = FALSE)
			}
		} else {
			if(pamFlag){
				progress$inc(0.01, detail = "Calculating MENTHUv2.0 scores for CRISPR sites...")
				
				# Function to calculate MENTHU v2.0 score
				menthuFrameFunc <- function(x){
					if(talFlag){
						progress$inc(1 / (nrow(pamSites) + nrow(talSites)))
					} else {
						progress$inc(1 /  nrow(pamSites))
					}
					
					return(calculateMenthu2(as.character(x), cutSite = 40, weight = 20, maxdbm = 5))
				}
				
				# Get the menthu scores
				menthuFrame <- as.data.frame(matrix(unlist(sapply(context, menthuFrameFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
				
				# Update progress bar
				progress$inc(0.01, detail = "Formatting CRISPR site results...")
				
				row.names(menthuFrame)  <- c()
				colnames(menthuFrame)   <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
				menthuFrame$menthuScore <- as.numeric(menthuFrame$menthuScore)
				
				# Get the critOne success count
				critOne <- nrow(menthuFrame[which(menthuFrame$frameShift != "NA"), ])
				
				# Get the critTwo success count
				critTwo <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5), ])
				
				# Get both success count
				critBoth <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5 & menthuFrame$frameShift != "NA"), ])
				
				pamSites <- unique(suppressMessages(plyr::join(pamSites, menthuFrame)))
				
				# Drop 0s
				# pamSites <- pamSites[which(pamSites$menthuScore > 0), ]
				
				# Format the output
				# Generate the 20bp CRISPR guide
				baseCrispr <- sapply(1:nrow(pamSites), 
														 function(x) if(pamSites$Orientation[x] == "forward"){
														 	substr(pamSites$seq[x], 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] - 19, 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x])     
														 } else {
														 	substr(reverseComplement(pamSites$seq[x]), 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] - 19, 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x])
														 })
				
				# Generate the PAM sequence
				pam        <- sapply(1:nrow(pamSites), 
														 function(x) (if(pamSites$Orientation[x] == "forward"){
														 	substr(pamSites$seq[x], 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + 1, 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))      
														 } else {
														 	substr(reverseComplement(pamSites$seq[x]), 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + 1, 
														 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))
														 }))
				
				# Get the CRISPR target sequence required to target this site, and bold the PAM
				crispr <- sapply(1:length(baseCrispr), function(x) paste0(baseCrispr[x], "<strong>", pam[x], "</strong>"))
				
				# Create data frame of the current results
				pamFormFrame  <- data.frame(Target_Sequence  = crispr, 
																		MENTHU_Score     = round(as.numeric(pamSites$menthuScore), digits = 2), 
																		Frame_Shift      = pamSites$frameShift,
																		Tool_Type        = pamSites$Target, 
																		Strand           = pamSites$Orientation, 
																		Exon_ID          = pamSites$Exon_Num, 
																		Cut_Location     = pamSites$CutIndex,
																		Top_Deletion     = pamSites$topDel,
																		stringsAsFactors = FALSE)
			}
			
			if(talFlag){
				# Update progress bar
				progress$inc(0.01, detail = "Calculating MENTHUv2.0 scores for TALEN sites...")
				
				menthuFrameTFunc <- function(x){
					# Increment progress
					if(pamFlag){
						progress$inc(1 / (nrow(talSites) + nrow(pamSites)))
						
					} else {
						progress$inc(1 / nrow(talSites))
						
					}
					
					# Return calculation
					return(calculateMenthu2(as.character(x), weight = 20, maxdbm = 5))
				}
				
				# Calculate slope competition on all the context
				menthuFrameT <- as.data.frame(matrix(unlist(sapply(contextT, menthuFrameTFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
				
				# Update progress bar
				progress$inc(0.01, detail = "Formatting TALEN site results...")
				
				# Clean up the resulting data frame
				colnames(menthuFrameT) <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
				rownames(menthuFrameT) <- c()
				
				# Merge slope frame and pamSites
				talSites <- unique(suppressMessages(plyr::join(talSites, menthuFrameT)))
				
				talenGenFunc <- function(talRow){
					dim  <- unlist(strsplit(talRow$Target, "/"))
					arm1 <- as.numeric(dim[1])
					spa1 <- as.numeric(dim[2])
					arm2 <- as.numeric(dim[3])
					
					armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
					spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
					armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
					
					return(paste0("<strong>", armL, "</strong>", spac, "<strong>", armR, "</strong>"))
				}
				
				# Generate the crispr target sequences
				talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
				
				# Clean the data frame
				row.names(talSites) <- c()
				
				# Create a new data frame of the results
				talFormFrame <- data.frame(Target_Sequence = talSites$talen,
																	 MENTHU_Score    = round(abs(as.numeric(talSites$menthuScore)), digits = 2),
																	 Frame_Shift     = talSites$frameShift,
																	 Tool_Type       = talSites$Target,
																	 Strand          = talSites$Orientation,
																	 Exon_ID         = talSites$Exon_Num,
																	 Cut_Location    = talSites$CutIndex,
																	 Top_Deletion    = talSites$topDel,
																	 stringsAsFactors = FALSE)
			}
		}
		
		#Return frame
		if(pamFlag && talFlag){
			return(list(rbind(pamFormFrame, talFormFrame), siteCount, siteCountC, critOne, critTwo, critBoth))
		} else if(pamFlag && !talFlag){
			return(list(pamFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
		} else {
			return(list(talFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
		}
		
	} else {
		return(1)
	}
}


#' calculateMENTHUGeneSeqGenBank
#'
#' @param pamList 
#' @param talenList 
#' @param gbFlag 
#' @param gbhFlag 
#' @param genbankInfo 
#' @param threshold 
#' @param firstExon 
#' @param exonTargetType 
#' @param exonStuff 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples

#calculateMENTHUGeneSeqGenBank <- function(pamList, cutDistList, wiggle = TRUE, wiggleRoom = 39, talenList, gbFlag, gbhFlag, talFlag,
#																					genbankInfo, threshold, firstExon, exonTargetType, exonStuff, progress, version){
#calculateMENTHUGeneSeqGenBank <- function(pamList, cutDistList, wiggle = TRUE, wiggleRoom = 39, talenList, gbFlag, gbhFlag, 
#																					genbankInfo, threshold, firstExon, exonTargetType, exonStuff, progress){
calculateMENTHUGeneSeqGenBank <- function(pamList, cutDistList, wiggle = TRUE, wiggleRoom = 39, talenList, gbFlag, gbhFlag, 
																					genbankInfo, firstExon, exonTargetType, exonStuff, progress){
	version <- 2
	require(plyr)
	
	# If NGG and NRG are both selected, only search for NGG to save time
	if("NGG" %in% pamList & "NRG" %in% pamList){
		pamList     <- pamList[-1]
		cutDistList <- cutDistList[-1]
	}
	
	# Update progress bar
	progress$inc(0.01, detail = "Processing GenBank accession...")
	
	# Get exon sequences and information through getExon
	exon <- getExon(genbankInfo, wiggle = TRUE, wiggleRoom = 39, gbFlag, exonTargetType, firstExon, exonStuff)
	
	# Get exon indices
	exonInfo <- exon[[1]]
	# Get the exon sequences
	exonSeq  <- exon[[2]]
	# Get the gene sequence
	geneSeq  <- exon[[3]]
	
	exonDF <- data.frame(Exon_Num         = exonInfo$exonNum,
											 exonStart        = exonInfo$start, 
											 exonEnd          = exonInfo$end, 
											 stringsAsFactors = FALSE)
	
	# If the user is using Cas:
	if(length(pamList) > 0){
		# Update progress bar
		progress$inc(0.01, detail = "Scanning for target sites...")
		
		if(length(exonInfo) > 0){
			# If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeq, 
													exonList   = exonInfo$exonNum, 
													exonStarts = exonInfo$start, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = TRUE, 
													wiggleRoom    = 39)
		} else {
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeqs,
													exonList   = "1",
													exonStarts = NULL, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wiggleRoom    = wiggleRoom)
		}
		
		siteCount <- nrow(pamSites)
		
		# Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
		# Update progress bar
		progress$inc(0.01, detail = "Pre-processing CRISPR target sites...")
		pamSites <- unique(suppressMessages(plyr::join(pamSites, exonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keep     <- sapply(1:nrow(pamSites), function(x) pamSites$CutIndex[x] %in% seq(from = pamSites$exonStart[x], to = pamSites$exonEnd[x], by = 1))
		pamSites <- pamSites[keep, ]
		
		# Identify sites with enough sequence context to do calculations
		pamSites$contextCondition[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		pamSites      <- pamSites[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		siteCountC <- nrow(pamSites)
		
		# Get the sequence context surrounding the cut site
		context <- unlist(lapply(1:nrow(pamSites), function (x) substr(geneSeq, pamSites$CutIndex[x] - 39, pamSites$CutIndex[x] + 40)))
		
		# Set the sequence context in the frame
		pamSites$seq = context
		
	} else {
		# If the user is NOT using Cas, set pamFlag to FALSE
		pamSites   <- 0
		pamFlag    <- FALSE
		
		siteCount  <- 0
		siteCountC <- 0
	}
	
	# Set a flag to be true if there are TALEN inputs
	talFlag <- talenList[1] != "" && talenList[2] != "" && talenList[3] != "" && talenList[4] != ""
	
	# If there are TALEN inputs
	if(talFlag){
		# Set the range flag to true
		rFlag <- TRUE
		
		# Set all exon starts to the exon starts in the input frame
		# Submit talen info to talPal
		# If there are exon inputs
		if(length(exonInfo > 0)){
			talSites <- talPal(exonSeq,
												 findCut    = TRUE,
												 wiggle     = TRUE,
												 wiggleRoom    = 39,
												 range      = rFlag, 
												 armin      = talenList[[1]], 
												 armax      = talenList[[2]], 
												 spamin     = talenList[[3]], 
												 spamax     = talenList[[4]], 
												 exonList   = exonInfo$exonNum,
												 exonStarts = exonInfo$start)
			
		} else {
			talSites <- talPal(exonSeq,
												 findCut    = TRUE,
												 wiggle     = TRUE,
												 wiggleRoom    = 39,
												 range      = rFlag,
												 armin      = talenList[[1]], 
												 armax      = talenList[[2]], 
												 spamin     = talenList[[3]], 
												 spamax     = talenList[[4]], 
												 exonStarts = NULL,
												 exonList   = "1")
		}
		
		# Update progress bar
		progress$inc(0.01, detail = "Pre-processing TALEN target sites...")
		talSites <- unique(suppressMessages(plyr::join(talSites, exonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keepT     <- sapply(1:nrow(talSites), function(x) talSites$CutIndex[x] %in% seq(from = talSites$exonStart[x], to = talSites$exonEnd[x], by = 1))
		talSites  <- talSites[keepT, ]
		
		# Identify sites with enough sequence context to do calculations
		talSites$contextCondition[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		talSites      <- talSites[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		# Get the sequence context surrounding the cut site
		contextT <- unlist(lapply(1:nrow(talSites), function (x) substr(geneSeq, talSites$CutIndex[x] - 39, talSites$CutIndex[x] + 40)))
		
		# Set the sequence context in the frame
		talSites$seq = contextT
		
	} else {
		# If TALENs are not used, set talSites list to empty
		talSites <- 0
	}
	
	# Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence  = as.character(), 
														MENTHU_Score     = as.numeric(), 
														Frame_Shift      = as.character(), 
														Tool_Type        = as.character(), 
														Strand           = as.character(), 
														Exon_ID          = as.numeric(), 
														Cut_Location     = as.integer(),
														stringsAsFactors = FALSE)
	if(version == 1){
		if(pamFlag){
			# Update progress bar
			progress$inc(0.01, detail = "Calculating MENTHUv1.0 scores for CRISPR sites...")
			
			slopeFrameFunc <- function(x){
				# Increment progress
				if(talFlag){
					progress$inc(1 / (nrow(pamSites) + nrow(talSites)))
					
				} else {
					progress$inc(1 / nrow(pamSites))
					
				}
				
				# Return calculation
				return(calculateSlopeCompetition(as.character(x), cutSite = 40, weight = 20, top = 10))
			}
			
			# Calculate slope competition on all the context
			slopeFrame <- as.data.frame(matrix(unlist(sapply(context, slopeFrameFunc)), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting CRISPR site results...")
			
			# Clean up the resulting data frame
			colnames(slopeFrame) <- c("seq", "microhomology_score", "OOF_Score", "slopeMH3Plus", "frameShift", "topDel")
			rownames(slopeFrame) <- c()
			
			# Merge slope frame and pamSites
			pamSites <- unique(suppressMessages(plyr::join(pamSites, slopeFrame)))
			
			# Generate the crispr target sequences
			pamSites$crispr <- sapply(1:nrow(pamSites), 
																function(x) substring((if(pamSites$Orientation[x] == "forward"){pamSites$seq[x]} else {reverseComplement(pamSites$seq[x])}),
																											40 - pamSites$CutDist[x] - 19,
																											40 - pamSites$CutDist[x] + nchar(pamSites$Target[x])))
			# Clean the data frame
			row.names(pamSites) <- c()
			
			# Create a new data frame of the results
			pamFormFrame <- data.frame(Target_Sequence = pamSites$crispr,
																 MENTHU_Score    = round(abs(as.numeric(pamSites$slopeMH3Plus)), digits = 2),
																 Frame_Shift     = pamSites$frameShift,
																 Tool_Type       = pamSites$Target,
																 Strand          = pamSites$Orientation,
																 Exon_ID         = pamSites$Exon_Num,
																 Cut_Location    = pamSites$CutIndex,
																 Top_Deletion    = pamSites$topDel,
																 stringsAsFactors = FALSE)
		}
		
		if(talFlag){
			# Update progress bar
			progress$inc(0.01, detail = "Calculating MENTHUv1.0 scores for TALEN sites...")
			
			slopeFrameTFunc <- function(x){
				# Increment progress
				if(pamFlag){
					progress$inc(1 / (nrow(talSites) + nrow(pamSites)))
					
				} else {
					progress$inc(1 / nrow(talSites))
					
				}
				
				# Return calculation
				return(calculateSlopeCompetition(as.character(x), weight = 20, top = 10))
			}
			
			# Calculate slope competition on all the context
			slopeFrameT <- as.data.frame(matrix(unlist(sapply(contextT, slopeFrameTFunc)), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting TALEN site results...")
			
			# Clean up the resulting data frame
			colnames(slopeFrameT) <- c("seq", "microhomology_score", "OOF_Score", "slopeMH3Plus", "frameShift", "topDel")
			rownames(slopeFrameT) <- c()
			
			# Merge slope frame and pamSites
			talSites <- unique(suppressMessages(plyr::join(talSites, slopeFrameT)))
			
			talenGenFunc <- function(talRow){
				dim  <- unlist(strsplit(talRow$Target, "/"))
				arm1 <- as.numeric(dim[1])
				spa1 <- as.numeric(dim[2])
				arm2 <- as.numeric(dim[3])
				
				armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
				spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
				armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
				
				return(paste0("<strong>", armL, "</strong>", spac, "<strong>", armR, "</strong>"))
			}
			
			# Generate the crispr target sequences
			talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
			
			# Clean the data frame
			row.names(talSites) <- c()
			
			# Create a new data frame of the results
			talFormFrame <- data.frame(Target_Sequence = talSites$talen,
																 MENTHU_Score    = round(abs(as.numeric(talSites$slopeMH3Plus)), digits = 2),
																 Frame_Shift     = talSites$frameShift,
																 Tool_Type       = talSites$Target,
																 Strand          = talSites$Orientation,
																 Exon_ID         = talSites$Exon_Num,
																 Cut_Location    = talSites$CutIndex,
																 Top_Deletion    = talSites$topDel,
																 stringsAsFactors = FALSE)
		}
	} else {
		if(pamFlag){
			progress$inc(0.01, detail = "Calculating MENTHUv2.0 scores for CRISPR sites...")
			
			# Function to calculate MENTHU v2.0 score
			menthuFrameFunc <- function(x){
				if(talFlag){
					progress$inc(1 / (nrow(pamSites) + nrow(talSites)))
				} else {
					progress$inc(1 /  nrow(pamSites))
				}
				
				return(calculateMenthu2(as.character(x), cutSite = 40, weight = 20, maxdbm = 5))
			}
			
			# Get the menthu scores
			menthuFrame <- as.data.frame(matrix(unlist(sapply(context, menthuFrameFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting CRISPR site results...")
			
			row.names(menthuFrame)  <- c()
			colnames(menthuFrame)   <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
			menthuFrame$menthuScore <- as.numeric(menthuFrame$menthuScore)
			
			# Get the critOne success count
			critOne <- nrow(menthuFrame[which(menthuFrame$frameShift != "NA"), ])
			
			# Get the critTwo success count
			critTwo <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5), ])
			
			# Get both success count
			critBoth <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5 & menthuFrame$frameShift != "NA"), ])
			
			pamSites <- unique(suppressMessages(plyr::join(pamSites, menthuFrame)))
			
			# Drop 0s
			pamSites <- pamSites[which(pamSites$menthuScore > 0), ]
			
			# Format the output
			# Generate the 20bp CRISPR guide
			baseCrispr <- sapply(1:nrow(pamSites), 
													 function(x) if(pamSites$Orientation[x] == "forward"){
													 	substr(pamSites$seq[x], 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] - 19, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x])     
													 } else {
													 	substr(reverseComplement(pamSites$seq[x]), 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] - 19, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x])
													 })
			
			# Generate the PAM sequence
			pam        <- sapply(1:nrow(pamSites), 
													 function(x) (if(pamSites$Orientation[x] == "forward"){
													 	substr(pamSites$seq[x], 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + 1, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))      
													 } else {
													 	substr(reverseComplement(pamSites$seq[x]), 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + 1, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))
													 }))
			
			# Get the CRISPR target sequence required to target this site, and bold the PAM
			crispr <- sapply(1:length(baseCrispr), function(x) paste0(baseCrispr[x], "<strong>", pam[x], "</strong>"))
			
			# Create data frame of the current results
			pamFormFrame  <- data.frame(Target_Sequence  = crispr, 
																	MENTHU_Score     = round(as.numeric(pamSites$menthuScore), digits = 2), 
																	Frame_Shift      = pamSites$frameShift,
																	Tool_Type        = pamSites$Target, 
																	Strand           = pamSites$Orientation, 
																	Exon_ID          = pamSites$Exon_Num, 
																	Cut_Location     = pamSites$CutIndex,
																	Top_Deletion     = pamSites$topDel,
																	stringsAsFactors = FALSE)
		}
		
		if(talFlag){
			# Update progress bar
			progress$inc(0.01, detail = "Calculating MENTHUv2.0 scores for TALEN sites...")
			
			menthuFrameTFunc <- function(x){
				# Increment progress
				if(pamFlag){
					progress$inc(1 / (nrow(talSites) + nrow(pamSites)))
					
				} else {
					progress$inc(1 / nrow(talSites))
					
				}
				
				# Return calculation
				return(calculateMenthu2(as.character(x), weight = 20, maxdbm = 5))
			}
			
			# Calculate slope competition on all the context
			menthuFrameT <- as.data.frame(matrix(unlist(sapply(contextT, menthuFrameTFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting TALEN site results...")
			
			# Clean up the resulting data frame
			colnames(menthuFrameT) <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
			rownames(menthuFrameT) <- c()
			
			# Merge slope frame and pamSites
			talSites <- unique(suppressMessages(plyr::join(talSites, menthuFrameT)))
			
			talenGenFunc <- function(talRow){
				dim  <- unlist(strsplit(talRow$Target, "/"))
				arm1 <- as.numeric(dim[1])
				spa1 <- as.numeric(dim[2])
				arm2 <- as.numeric(dim[3])
				
				armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
				spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
				armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
				
				return(paste0("<strong>", armL, "</strong>", spac, "<strong>", armR, "</strong>"))
			}
			
			# Generate the crispr target sequences
			talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
			
			# Clean the data frame
			row.names(talSites) <- c()
			
			# Create a new data frame of the results
			talFormFrame <- data.frame(Target_Sequence = talSites$talen,
																 MENTHU_Score    = round(abs(as.numeric(talSites$menthuScore)), digits = 2),
																 Frame_Shift     = talSites$frameShift,
																 Tool_Type       = talSites$Target,
																 Strand          = talSites$Orientation,
																 Exon_ID         = talSites$Exon_Num,
																 Cut_Location    = talSites$CutIndex,
																 Top_Deletion    = talSites$topDel,
																 stringsAsFactors = FALSE)
		}
	}
	
	# Return frame
	if(pamFlag && talFlag){
		return(list(rbind(pamFormFrame, talFormFrame), siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(pamFlag && !talFlag){
		return(list(pamFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	} else {
		return(list(talFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	}
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> parent of f2dd323... Revert "Added functions for Ensembl support; fixed issues with generating 20 bp gRNA base; fixed several parentheses"
=======
>>>>>>> parent of f2dd323... Revert "Added functions for Ensembl support; fixed issues with generating 20 bp gRNA base; fixed several parentheses"

#' calculateMENTHUEnsembl
#'
#' @param pamList 
#' @param talenList 
#' @param gbFlag 
#' @param gbhFlag 
#' @param genbankInfo 
#' @param threshold 
#' @param firstExon 
#' @param exonTargetType 
#' @param exonStuff 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples

#stuff    <- calculateMENTHUEnsembl(pams, cutDistances, wiggle = TRUE, wiggleRoom = 39, talenList, 
#																	 info, input$firstExon, input$exonTargetType, exonStuff, progress)

#calculateMENTHUEnsembl <- function(pamList, cutDistList, wiggle = TRUE, wiggleRoom = 39, talenList, talFlag,
#																					genbankInfo, threshold, firstExon, exonTargetType, exonStuff, progress, version){
#calculateMENTHUEnsembl <- function(pamList, cutDistList, wiggle = TRUE, wiggleRoom = 39, talenList, 
#																					genbankInfo, threshold, firstExon, exonTargetType, exonStuff, progress){
calculateMENTHUEnsembl <- function(pamList, cutDistList, wiggle = TRUE, wiggleRoom = 39, talenList, ensemblInfo, exonStuff, progress){
	version <- 2
	require(plyr)
	
	# If NGG and NRG are both selected, only search for NGG to save time
	if("NGG" %in% pamList & "NRG" %in% pamList){
		pamList     <- pamList[-1]
		cutDistList <- cutDistList[-1]
	}
	
	# Update progress bar
	progress$inc(0.01, detail = "Processing Ensembl sites...")
	
<<<<<<< HEAD
<<<<<<< HEAD
	# Generate a subset of exons
=======
>>>>>>> parent of f2dd323... Revert "Added functions for Ensembl support; fixed issues with generating 20 bp gRNA base; fixed several parentheses"
=======
>>>>>>> parent of f2dd323... Revert "Added functions for Ensembl support; fixed issues with generating 20 bp gRNA base; fixed several parentheses"
	exonSubset <- ensemblInfo[which(as.numeric(ensemblInfo$rank) %in% as.numeric(exonStuff)), ]
	
	# Get the exon sequences
	exonSeq  <- exonSubset$sequence
	
	# Create a data frame to hold information about the exon sequence
	exonDF <- data.frame(Exon_Num         = as.numeric(exonSubset$rank),
											 absStart         = as.numeric(exonSubset$contextStart),
											 absEnd           = as.numeric(exonSubset$contextEnd),
											 exonStart        = 0 + as.numeric(exonSubset$exp5),
											 exonEnd          = nchar(exonSubset$sequence) - as.numeric(exonSubset$exp3),
											 seq              = exonSeq,
											 stringsAsFactors = FALSE)
	
	# If the user is using Cas:
	if(length(pamList) > 0){
		# Update progress bar
		progress$inc(0.01, detail = "Scanning for target sites...")
		
		#if(length(exonStuff) > 0){
			# If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeq, 
													exonList   = exonDF$Exon_Num, 
													exonStarts = NULL, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wiggleRoom    = wiggleRoom)
		#} else {
		#	pamSites <- pamScan(pamList, 
		#											cutDistList, 
		#											exonSeqs,
		#											exonList   = "1",
		#											exonStarts = NULL, 
		#											findCut    = TRUE, 
		#											type       = "cas9", 
		#											wiggle     = wiggle, 
		#											wiggleRoom    = wiggleRoom)
		#}
		
		siteCount <- nrow(pamSites)
		
		# Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
		# Update progress bar
		progress$inc(0.01, detail = "Pre-processing CRISPR target sites...")
		pamSites        <- unique(suppressMessages(plyr::join(pamSites, exonDF, by = 'Exon_Num')))
		
		# THIS IS WRONG
		pamSites$absCut <- pamSites$CutIndex + pamSites$absStart
		
		# Drop target sites where the cut site is not within the exon boundaries
		keep     <- sapply(1:nrow(pamSites), function(x) pamSites$CutIndex[x] %in% seq(from = pamSites$exonStart[x], to = pamSites$exonEnd[x], by = 1))
		pamSites <- pamSites[keep, ]
		
		# Identify sites with enough sequence context to do calculations
		pamSites$contextCondition[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(pamSites$seq)))  ] <- TRUE
		pamSites      <- pamSites[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(pamSites$seq))), ]
		
		siteCountC <- nrow(pamSites)
		
		# Get the sequence context surrounding the cut site
		context <- unlist(lapply(1:nrow(pamSites), function (x) substr(pamSites$seq[x], pamSites$CutIndex[x] - 39, pamSites$CutIndex[x] + 40)))
		
		# Set the sequence context in the frame
		pamSites$seq = context
		
	} else {
		# If the user is NOT using Cas, set pamFlag to FALSE
		pamSites   <- 0
		pamFlag    <- FALSE
		
		siteCount  <- 0
		siteCountC <- 0
	}
	
	# Set a flag to be true if there are TALEN inputs
	talFlag <- talenList[1] != "" && talenList[2] != "" && talenList[3] != "" && talenList[4] != ""
	
	# If there are TALEN inputs
	if(talFlag){
		# Set the range flag to true
		rFlag <- TRUE
		
		# Set all exon starts to the exon starts in the input frame
		# Submit talen info to talPal
		# If there are exon inputs
		if(length(exonInfo) > 0){
			talSites <- talPal(exonSeq,
												 findCut    = TRUE,
												 wiggle     = TRUE,
												 wiggleRoom    = 39,
												 range      = rFlag, 
												 armin      = talenList[[1]], 
												 armax      = talenList[[2]], 
												 spamin     = talenList[[3]], 
												 spamax     = talenList[[4]], 
												 exonList   = exonInfo$exonNum,
												 exonStarts = exonInfo$start)
			
		} else {
			talSites <- talPal(exonSeq,
												 findCut    = TRUE,
												 wiggle     = TRUE,
												 wiggleRoom    = 39,
												 range      = rFlag,
												 armin      = talenList[[1]], 
												 armax      = talenList[[2]], 
												 spamin     = talenList[[3]], 
												 spamax     = talenList[[4]], 
												 exonStarts = NULL,
												 exonList   = "1")
		}
		
		# Update progress bar
		progress$inc(0.01, detail = "Pre-processing TALEN target sites...")
		talSites <- unique(suppressMessages(plyr::join(talSites, exonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keepT     <- sapply(1:nrow(talSites), function(x) talSites$CutIndex[x] %in% seq(from = talSites$exonStart[x], to = talSites$exonEnd[x], by = 1))
		talSites  <- talSites[keepT, ]
		
		# Identify sites with enough sequence context to do calculations
		talSites$contextCondition[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		talSites      <- talSites[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		# Get the sequence context surrounding the cut site
		contextT <- unlist(lapply(1:nrow(talSites), function (x) substr(geneSeq, talSites$CutIndex[x] - 39, talSites$CutIndex[x] + 40)))
		
		# Set the sequence context in the frame
		talSites$seq = contextT
		
	} else {
		# If TALENs are not used, set talSites list to empty
		talSites <- 0
	}
	
	# Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence  = as.character(), 
														MENTHU_Score     = as.numeric(), 
														Frame_Shift      = as.character(), 
														Tool_Type        = as.character(), 
														Strand           = as.character(), 
														Exon_ID          = as.numeric(), 
														Cut_Location     = as.integer(),
														stringsAsFactors = FALSE)
	if(version == 1){
		if(pamFlag){
			# Update progress bar
			progress$inc(0.01, detail = "Calculating MENTHUv1.0 scores for CRISPR sites...")
			
			slopeFrameFunc <- function(x){
				# Increment progress
				if(talFlag){
					progress$inc(1 / (nrow(pamSites) + nrow(talSites)))
					
				} else {
					progress$inc(1 / nrow(pamSites))
					
				}
				
				# Return calculation
				return(calculateSlopeCompetition(as.character(x), cutSite = 40, weight = 20, top = 10))
			}
			
			# Calculate slope competition on all the context
			slopeFrame <- as.data.frame(matrix(unlist(sapply(context, slopeFrameFunc)), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting CRISPR site results...")
			
			# Clean up the resulting data frame
			colnames(slopeFrame) <- c("seq", "microhomology_score", "OOF_Score", "slopeMH3Plus", "frameShift", "topDel")
			rownames(slopeFrame) <- c()
			
			# Merge slope frame and pamSites
			pamSites <- unique(suppressMessages(plyr::join(pamSites, slopeFrame)))
			
			# Generate the crispr target sequences
			pamSites$crispr <- sapply(1:nrow(pamSites), 
																function(x) substring((if(pamSites$Orientation[x] == "forward"){pamSites$seq[x]} else {reverseComplement(pamSites$seq[x])}),
																											40 - pamSites$CutDist[x] - 19,
																											40 - pamSites$CutDist[x] + nchar(pamSites$Target[x])))
			# Clean the data frame
			row.names(pamSites) <- c()
			
			# Create a new data frame of the results
			pamFormFrame <- data.frame(Target_Sequence = pamSites$crispr,
																 MENTHU_Score    = round(abs(as.numeric(pamSites$slopeMH3Plus)), digits = 2),
																 Frame_Shift     = pamSites$frameShift,
																 Tool_Type       = pamSites$Target,
																 Strand          = pamSites$Orientation,
																 Exon_ID         = pamSites$Exon_Num,
																 Cut_Location    = pamSites$CutIndex,
																 Top_Deletion    = pamSites$topDel,
																 stringsAsFactors = FALSE)
		}
		
		if(talFlag){
			# Update progress bar
			progress$inc(0.01, detail = "Calculating MENTHUv1.0 scores for TALEN sites...")
			
			slopeFrameTFunc <- function(x){
				# Increment progress
				if(pamFlag){
					progress$inc(1 / (nrow(talSites) + nrow(pamSites)))
					
				} else {
					progress$inc(1 / nrow(talSites))
					
				}
				
				# Return calculation
				return(calculateSlopeCompetition(as.character(x), weight = 20, top = 10))
			}
			
			# Calculate slope competition on all the context
			slopeFrameT <- as.data.frame(matrix(unlist(sapply(contextT, slopeFrameTFunc)), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting TALEN site results...")
			
			# Clean up the resulting data frame
			colnames(slopeFrameT) <- c("seq", "microhomology_score", "OOF_Score", "slopeMH3Plus", "frameShift", "topDel")
			rownames(slopeFrameT) <- c()
			
			# Merge slope frame and pamSites
			talSites <- unique(suppressMessages(plyr::join(talSites, slopeFrameT)))
			
			talenGenFunc <- function(talRow){
				dim  <- unlist(strsplit(talRow$Target, "/"))
				arm1 <- as.numeric(dim[1])
				spa1 <- as.numeric(dim[2])
				arm2 <- as.numeric(dim[3])
				
				armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
				spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
				armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
				
				return(paste0("<strong>", armL, "</strong>", spac, "<strong>", armR, "</strong>"))
			}
			
			# Generate the crispr target sequences
			talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
			
			# Clean the data frame
			row.names(talSites) <- c()
			
			# Create a new data frame of the results
			talFormFrame <- data.frame(Target_Sequence = talSites$talen,
																 MENTHU_Score    = round(abs(as.numeric(talSites$slopeMH3Plus)), digits = 2),
																 Frame_Shift     = talSites$frameShift,
																 Tool_Type       = talSites$Target,
																 Strand          = talSites$Orientation,
																 Exon_ID         = talSites$Exon_Num,
																 Cut_Location    = talSites$CutIndex,
																 Top_Deletion    = talSites$topDel,
																 stringsAsFactors = FALSE)
		}
	} else {
		if(pamFlag){
			progress$inc(0.01, detail = "Calculating MENTHUv2.0 scores for CRISPR sites...")
			
			# Function to calculate MENTHU v2.0 score
			menthuFrameFunc <- function(x){
				if(talFlag){
					progress$inc(1 / (nrow(pamSites) + nrow(talSites)))
				} else {
					progress$inc(1 /  nrow(pamSites))
				}
				
				return(calculateMenthu2(as.character(x), cutSite = 40, weight = 20, maxdbm = 5))
			}
			
			# Get the menthu scores
			menthuFrame <- as.data.frame(matrix(unlist(sapply(context, menthuFrameFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting CRISPR site results...")
			
			row.names(menthuFrame)  <- c()
			colnames(menthuFrame)   <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
			menthuFrame$menthuScore <- as.numeric(menthuFrame$menthuScore)
			
			# Get the critOne success count
			critOne <- nrow(menthuFrame[which(menthuFrame$frameShift != "NA"), ])
			
			# Get the critTwo success count
			critTwo <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5), ])
			
			# Get both success count
			critBoth <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5 & menthuFrame$frameShift != "NA"), ])
			
			pamSites <- unique(suppressMessages(plyr::join(pamSites, menthuFrame)))
			
			# Drop 0s
			pamSites <- pamSites[which(pamSites$menthuScore > 0), ]
			
			# Format the output
			# Generate the 20bp CRISPR guide
			baseCrispr <- sapply(1:nrow(pamSites), 
													 function(x) if(pamSites$Orientation[x] == "forward"){
													 	substr(pamSites$seq[x], 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] - 19, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x])     
													 } else {
													 	substr(reverseComplement(pamSites$seq[x]), 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] - 19, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x])
													 })
			
			# Generate the PAM sequence
			pam        <- sapply(1:nrow(pamSites), 
													 function(x) (if(pamSites$Orientation[x] == "forward"){
													 	substr(pamSites$seq[x], 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + 1, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))      
													 } else {
													 	substr(reverseComplement(pamSites$seq[x]), 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + 1, 
													 				 (nchar(pamSites$seq[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))
													 }))
			
			# Get the CRISPR target sequence required to target this site, and bold the PAM
			crispr <- sapply(1:length(baseCrispr), function(x) paste0(baseCrispr[x], "<strong>", pam[x], "</strong>"))
			
			# Create data frame of the current results
			pamFormFrame  <- data.frame(Target_Sequence  = crispr, 
																	MENTHU_Score     = round(as.numeric(pamSites$menthuScore), digits = 2), 
																	Frame_Shift      = pamSites$frameShift,
																	Tool_Type        = pamSites$Target, 
																	Strand           = pamSites$Orientation, 
																	Exon_ID          = pamSites$Exon_Num, 
																	Cut_Location     = pamSites$CutIndex,
																	Top_Deletion     = pamSites$topDel,
																	stringsAsFactors = FALSE)
		}
		
		if(talFlag){
			# Update progress bar
			progress$inc(0.01, detail = "Calculating MENTHUv2.0 scores for TALEN sites...")
			
			menthuFrameTFunc <- function(x){
				# Increment progress
				if(pamFlag){
					progress$inc(1 / (nrow(talSites) + nrow(pamSites)))
					
				} else {
					progress$inc(1 / nrow(talSites))
					
				}
				
				# Return calculation
				return(calculateMenthu2(as.character(x), weight = 20, maxdbm = 5))
			}
			
			# Calculate slope competition on all the context
			menthuFrameT <- as.data.frame(matrix(unlist(sapply(contextT, menthuFrameTFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
			
			# Update progress bar
			progress$inc(0.01, detail = "Formatting TALEN site results...")
			
			# Clean up the resulting data frame
			colnames(menthuFrameT) <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
			rownames(menthuFrameT) <- c()
			
			# Merge slope frame and pamSites
			talSites <- unique(suppressMessages(plyr::join(talSites, menthuFrameT)))
			
			talenGenFunc <- function(talRow){
				dim  <- unlist(strsplit(talRow$Target, "/"))
				arm1 <- as.numeric(dim[1])
				spa1 <- as.numeric(dim[2])
				arm2 <- as.numeric(dim[3])
				
				armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
				spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
				armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
				
				return(paste0("<strong>", armL, "</strong>", spac, "<strong>", armR, "</strong>"))
			}
			
			# Generate the crispr target sequences
			talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
			
			# Clean the data frame
			row.names(talSites) <- c()
			
			# Create a new data frame of the results
			talFormFrame <- data.frame(Target_Sequence = talSites$talen,
																 MENTHU_Score    = round(abs(as.numeric(talSites$menthuScore)), digits = 2),
																 Frame_Shift     = talSites$frameShift,
																 Tool_Type       = talSites$Target,
																 Strand          = talSites$Orientation,
																 Exon_ID         = talSites$Exon_Num,
																 Cut_Location    = talSites$CutIndex,
																 Top_Deletion    = talSites$topDel,
																 stringsAsFactors = FALSE)
		}
	}
	
	# Return frame
	if(pamFlag && talFlag){
		return(list(rbind(pamFormFrame, talFormFrame), siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(pamFlag && !talFlag){
		return(list(pamFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	} else {
		return(list(talFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	}
}


<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> parent of 6a4cadd... Revert "Added comment for clarity"
=======
>>>>>>> parent of f2dd323... Revert "Added functions for Ensembl support; fixed issues with generating 20 bp gRNA base; fixed several parentheses"
=======
>>>>>>> parent of f2dd323... Revert "Added functions for Ensembl support; fixed issues with generating 20 bp gRNA base; fixed several parentheses"
#' convertToNumeric
#'
#' This function takes an input string of comma- or whitespace-separated numbers and converts to a numeric vector
#' @param characterStringOfNumbers 
#'
#' @return
#' @export
#'
#' @examples

convertToNumeric <- function(characterStringOfNumbers){
	# Split the string by ','
	splittedString <- unlist(strsplit(characterStringOfNumbers, "[\\, |\\,| ]+"))
	
	# If there is more than one object in the list
	if(length(splittedString) > 1){
		# Find instances where a range was specified
		rangeLocs <- grep("-", splittedString, fixed = TRUE)
		
		# If more than one range
		if(length(rangeLocs) > 1){
			rangeEnds <- strsplit(splittedString[rangeLocs], "-")
			numberList <- splittedString[-rangeLocs]
			
			# Generate a sequence in the range
			for(i in 1:length(rangeEnds)){
				numberList <- c(numberList, seq(from = rangeEnds[[i]][1], to = rangeEnds[[i]][2]))
			}
		} else {
			numberList <- splittedString
		}
		
		numberList <- as.numeric(unlist(numberList))
		numberList <- sort(numberList)
		
	} else {
		numberList <- as.numeric(splittedString)
	}
	
	return(numberList)
}


#' distStitch
#'
#' This function matches up cut indices with PAMs; e.g., all pre-computed PAMs currently cut 3bp
#' upstream of the PAM, so a '-3' is generated for each pre-computed PAM. This is then added on to
#' any custom PAM sites, with custom cut distances
#' 
#' @param pamList 
#' @param custList 
#'
#' @return
#' @export
#'
#' @examples
#' 

distStitch <- function(pamList, custList){
	# If there are pre-gen PAMs used...
	if(pamList != ""){
		# Repeat -3 for each pre-computed PAM (the ones currently on the list all cut 3 bp upstream of the PAM)
		pamDistList  <- rep(-3, length(pamList))
		# Split the custom cut distance list
		custDistList <- strsplit(custList, "[\\, |\\,| ]+")
		# Add the cut distances to a list
		distList     <- unlist(c(pamDistList, custDistList))
		
	} else {
		# Add the cut distances to a list (after splitting them)
		distList     <- unlist(strsplit(custList, "[\\, |\\,| ]+"))
	}
	
	return(as.numeric(distList))
}


#' pamStitch 
#' 
#' This function adds custom PAMs to pre-generated PAMs
#'
#' @param pamList 
#' @param custPamList 
#'
#' @return
#' @export
#'
#' @examples

pamStitch <- function(pamList, custPamList){
	# If there are pre-gen PAMs used...
	if(pamList != ""){
		# Split the custom PAM list
		custList   <- strsplit(custPamList, "[\\, |\\,| ]+")
		# Add to the pre-gen PAM list
		customList <- unlist(c(pamList, custList))
	} else {
		# Split the custom PAM list
		customList <- unlist(strsplit(custPamList, "[\\, |\\,| ]+"))
	}
	
	return(customList)
}


####Exon Handler for Custom Exon Input####
exonHandler <- function(exonRHandsonTable){
	exonTable <- hot_to_r(exonRHandsonTable)
	exonDF <- exonTable[apply(exonTable, MARGIN = 1, function(x) any(x > 0)),]
	
	return(exonDF)
}


#' window
#'  
#' This function gets the necessary sequence context to do MENTHU calculations
#' 
#' @examples 
#' 
#' @export 
#' 

window <- function(sequence, position, winSize = 80) {
	
	if(position < (winSize / 2)){
		return("")
		
	} else {
		return(substring(sequence, ((position - ((winSize + (winSize %% 2)) / 2) - 1)), (position + ((winSize + (winSize %% 2)) / 2))))
	}
	
}


#' reverse
#'
#' This function takes a string as input and reverses the order of the string
#'
#' @param seq A string to reverse
#'
#' @return revSeq The seq string in reverse
#'
#' @examples
#' reverse("123456")
#'
#'
#' @export

reverse <- function(seq){
	UseMethod("reverse", seq)
}

reverse.default <- function(seq){
	stop("Error: Cannot reverse objects that are not character strings or integers. Please check input sequence.")
}


reverse.character <- function(seq){
	revSeq <- seq
	
	for(i in 1:nchar(seq)){
		curL <- substr(seq, i, i)
		substr(revSeq, (nchar(seq) + 1 - i), (nchar(seq) + 1 - i)) <- curL
	}
	
	return(revSeq)
}

reverse.numeric <- function(seq){
	charSeq <- as.character(seq)
	revSeq <- charSeq
	revSeq <- as.numeric(reverse.character(charSeq))
	return(revSeq)
}

#' complement
#'
#' This function takes a DNA or RNA sequence as input (along with a parameter specifying the type of sequence) 
#' and outputs the complement of the input sequence. E.g., "ATTG" will return "TAAC" if type = "DNA" and "UAAC" if type = "RNA"
#'
#' @param seq A DNA or RNA sequence from which to generate a complement string
#' @param type (Now defunct; currently, all IUPAC 1-letter nucleotide codes are supported.) Default is "DNA"; a DNA sequence can only contain 
#' "A", "C", "G", or "T" for the purposes of complement(). The other option is "RNA"; an RNA sequence can only contain 
#' "A", "C", "G", or "U" for the purposes of complement().
#'
#' @return compSeq The complement of the input sequence
#' @export
#'
#' @examples
#' seq <- "AAAATGGCGAAG"
#' type <- "DNA"
#' complement(seq, type)

complement <- function(seq, type){
	UseMethod("complement", seq)
}

complement.default <- function(seq, type){
	#Prevent attempts to complement objects that are not sequences
	stop("Error: Cannot complement a non-character vector object. Please check input sequence.")
}

complement.character <- function(seq, type){
	compSeq <- seq
	
	fromVal <- c("A", "C", "G", "T", "a", "c", "g", "t", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v", "n")
	toVal   <- c("T", "G", "C", "A", "t", "g", "c", "a", "Y", "R", "S", "W", "M", "K", "V", "H", "D", "B", "N", "y", "r", "s", "w", "m", "k", "v", "h", "d", "b", "n")
	
	
	compSeq <- plyr::mapvalues(unlist(strsplit(compSeq, split = "")),
														 from         = fromVal,
														 to           = toVal,
														 warn_missing = FALSE)
	compSeq <- paste(compSeq, collapse = "")
	return(compSeq)
}

complement.list <- function(seq, type){
	retList <- list()
	for(i in seq){
		retList <- c(retList, complement(i, type))
	}
	return(retList)
}

#' reverseComplement
#'
#' This function takes a DNA or RNA sequence as input and outputs the reverse complement of the sequence.
#'
#' @param seq A character vector from which to generate a reverse complement.
#' @param type This is deprecated, as 'complement' now deals with all IUPAC 1-letter codes. 
#' Default is "DNA"; allowed characters are "A", "C", "G", and "T" (case insensitive). 
#' Other option is "RNA"; allowed characters are "A", "C", "G", and "U" (case insensitive.)
#'
#' @return seqRevComp The reverse complement of the input sequence
#' @export
#'
#' @examples
#' dnaSeq <- "AATGCC"
#' reverseComplement(dnaSeq)
#' rnaSeq <- "UUAGCC"
#' reverseComplement(rnaSeq, type = "RNA")

reverseComplement <- function(seq, type = "DNA"){
	UseMethod("reverseComplement", seq)
}

reverseComplement.default <- function(seq, type = "DNA"){
	stop("Error: Input sequence is not a character vector. Please check input sequence.")
}

reverseComplement.character <- function(seq, type = "DNA"){
	#Reverse the sequence
	seqRev <- reverse(seq)
	#Get the complement of the reversed sequence
	seqRevComp <- complement(seqRev, type = "DNA")
	return(seqRevComp)
}

reverseComplement.list <- function(seq, type = "DNA"){
	retList <- list()
	for(i in seq){
		retList <- c(retList, reverseComplement.character(i))
	}
	return(unlist(retList))
}
