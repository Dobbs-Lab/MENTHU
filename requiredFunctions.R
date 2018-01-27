####Required Function#######

#' Title
#'
#' @param casList 
#' @param wiggle 
#' @param wigRoom 
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

calculateMENTHUGeneSeq <- function(casList, cutDistList, wiggle = TRUE, wigRoom = 39, geneSeq, threshold, exonDF, progress, armin, armax, spamin, spamax){
	require(Biostrings)
	
	#Rename the input casList to pamList
	pamList <- casList
	
	#Create an empty list to hold exon sequences
	set <- NULL
	exonSeqs <- NULL
	
	#Subset geneSeq to exon sequences
	#Deal with the case in which exons are specified and extra context should be included to examine sites where the gRNA would run off the exon
	if((wiggle == TRUE) && (class(exonDF) == "data.frame")){
		for(i in 1:nrow(exonDF)){
			#Ensure extra context doesn't run off the beginning of the whole sequence
			if((exonDF$exonStart[i] - wigRoom) < 1){
				exStart <- 1
			} else {
				exStart <- exonDF$exonStart[i] - wigRoom
			}
			
			#Ensure extra context doesn't run off the end of the whole sequence
			if(exonDF$exonEnd[i] + wigRoom > nchar(geneSeq)){
				exEnd <- nchar(geneSeq)
			} else {
				exEnd <- exonDF$exonEnd[i] + wigRoom
			}
			
			#Create new exonSeqs set with 'fixed' context
			set <- c(set, DNAString(geneSeq)[exStart:exEnd])
			exonSeqs <- DNAStringSet(set)
		}
		
		#Deal with the case in which exons are specified but extra context should not be included
	} else if((!wiggle) && (class(exonDF) == "data.frame")){
		exonSeqs <- DNAStringSet(substring(rep(as.character(geneSeq), nrow(exonDF)), exonDF$exonStart, exonDF$exonEnd))
		
	} else {
		#In the case where exons are not specified, treat the whole gene sequence as a single exon
		exonSeqs <- DNAStringSet(geneSeq)
	}
	
	#Update progress bar
	progress$inc(0.01, detail = "Scanning for target sites...")
	
	#If the user is using Cas:
	if(length(pamList) > 0){
		#If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
		if(class(exonDF) == "data.frame"){
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeqs, 
													exonList   = row.names(exonDF), 
													exonStarts = exonDF$exonStart, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wigRoom    = wigRoom)
		} else {
			pamSites <- pamScan(pamList, 
													cutDistList, 
													exonSeqs, 
													exonStarts = NULL, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wigRoom    = wigRoom,
													exonList   = "1")
		}
		
		#Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
	} else {
		
		#If the user is NOT using CRISPR/Cas system, set pamFlag to FALSE
		pamSites <- 0
		pamFlag <- FALSE
	}
	
	#Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence = as.character(), 
														MENTHU_Score = as.numeric(), 
														Frame_Shift = as.character(), 
														Tool_Type = as.character(), 
														Strand = as.character(), 
														Exon_ID = as.numeric(), 
														Cut_Location = as.integer())
	
	#Set a flag to be true if there are TALEN inputs
	talFlag <- armin != "" & armax != "" & spamin != "" & spamax != ""
	
	#If there are TALEN inputs
	if(talFlag){
		#Set the range flag to true
		rFlag <- TRUE
		
		#If there are exon inputs
		if(class(exonDF) == "data.frame"){
			#Set all exon starts to the exon starts in the input frame
			#Submit talen info to talPal
			talSites <- talPal(exonSeqs,
												 findcut = TRUE,
												 wiggle = TRUE,
												 wigRoom = 39,
												 range      = rFlag, 
												 armin      = armin, 
												 armax      = armax, 
												 spamin     = spamin, 
												 spamax     = spamax, 
												 exonStarts = exonDF$exonStart,
												 exonList   = row.names(exonDF))
		} else {
			#If there are no exon inputs, make exon start null
			#Submit talen info to talPal
			talSites <- talPal(exonSeqs,
												 findcut = TRUE,
												 wiggle = TRUE,
												 wigRoom = 39,
												 range      = rFlag, 
												 armin      = armin, 
												 armax      = armax, 
												 spamin     = spamin, 
												 spamax     = spamax, 
												 exonStarts = NULL,
												 exonList   = "1")
		}
		
	} else {
		#If TALENs are not used, set talSites list to empty
		talSites <- ""
	}
	
	#If there is no exon data input, create a dummy data frame to have the 'exon' be the entire gene sequence, starting on nt 1
	if(class(exonDF) != "data.frame"){
		exonDF <- data.frame(exonStart = 1, exonEnd = nchar(geneSeq), stringsAsFactors = FALSE)
	}
	
	#If the user is using PAMS:
	if(pamFlag){
		
		#PAM LEVEL
		for(i in 1:length(pamSites[])){
			toolTypeI <- names(pamSites)[i]
			
			#STRAND LEVEL
			for(j in 1:length(pamSites[[i]][])){
				
				#EXON LEVEL
				for(k in 1:length(pamSites[[i]][[j]][])){
					#Get the number of the exon
					exonNum <- as.numeric(gsub("Exon ", "", names(pamSites[[i]][[j]])[k]))
					
					#SITE LEVEL
					for(l in 1:length(pamSites[[i]][[j]][[k]][])){
						
						#Make the process bar increment for every potential Cas site
						if(talFlag == FALSE){
							progress$inc(1/(length(unlist(pamSites))), detail = paste("Processing ", names(pamSites[i])))
						} else {
							progress$inc(1/(length(unlist(pamSites)) + sum(lengths(talSites[[1]][]))), detail = paste("Processing ", names(pamSites[i])))
						}
						
						if(!is.na(pamSites[[i]][[j]][[k]][[l]])){
							inExon <- FALSE
							exonId <- 0
							
							#Ignore target sites where the cut site is not within the current exon
							for(m in 1:nrow(exonDF)){
								if((pamSites[[i]][[j]][[k]][[l]] >= exonDF$exonStart[m]) && 
									 (pamSites[[i]][[j]][[k]][[l]] <= exonDF$exonEnd[m])){
									inExon <- TRUE
									exonId <- m
								}
							}
							
							if(exonId == exonNum){
								#If on the positive strand...
								if(j == 1){
									strandId <- "forward"
								} else if (j == 2){
									#If on the negative strand...
									strandId <- "complement" 
								}
								
								#Boolean value determining if there is enough sequence context (40 bp upstream, 40 bp downstream) of the cut site
								#in order to run deletion calculation
								contextCondition <- ((pamSites[[i]][[j]][[k]][[l]] > 40) && 
																		 	(pamSites[[i]][[j]][[k]][[l]] + 40) < nchar(geneSeq))
								
								#Check if there's enough context...
								if(!is.na(pamSites[[i]][[j]][[k]][[l]]) && !is.null(pamSites[[i]][[j]][[k]][[l]])){
									
									if(contextCondition){
										#Get the sequence context surrounding the cut site
										context <- window(DNAString(geneSeq), pamSites[[i]][[j]][[k]][[l]], 80)
										
										#Calculate competition for the sequence generated by window()
										slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
										
										#If the score is greater than the threshold, report it
										if(abs(slopeFrame$slopeMH3Plus) >= threshold){
											
											#Get the CRISPR target sequence required to target this site
											if(strandId == "forward"){
												crispr <- substr(slopeFrame$seq, 
																				 (40 + -cutDistList[i] - 19), 
																				 (40 + -cutDistList[i] + nchar(toolTypeI)))
											} else if(strandId == "complement"){
												crispr <- as.character(substr(Biostrings::reverseComplement(DNAString(slopeFrame$seq)), 
																											(40 + -cutDistList[i] - 19), 
																											(40 + -cutDistList[i] + nchar(toolTypeI))))
											}
											
											#Create data frame of the current results
											formFrame  <- data.frame(Target_Sequence = crispr, 
																							 MENTHU_Score    = round(abs(slopeFrame$slopeMH3Plus), digits = 2), 
																							 Frame_Shift     = slopeFrame$frameShift,
																							 Tool_Type       = toolTypeI, 
																							 Strand          = strandId, 
																							 Exon_ID         = exonNum, 
																							 Cut_Location    = as.integer(pamSites[[i]][[j]][[k]][[l]], digits = 0))
											
											#Bind the current results to the running list
											menthuFrame <- rbind(menthuFrame, formFrame)
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	if(talFlag){
		
		#For TALENs
		toolTypeI <- "TALEN"
		
		for(o in 1:length(talSites[[1]][])){
			
			if(length(talSites[[1]][[o]] > 0)){
				#Get the number of the exon
				exonNum <- names(talSites[[2]][o])
				exonNumAct <- as.numeric(gsub("Exon ", "", exonNum))
				
				for(p in 1:length(talSites[[1]][[o]])){
					#Only do the calculation of the TALEN cut site if it is in the current exon
					if(length(talSites[[1]][[o]][[p]][]) > 0){
						#Progress bar for going through sites
						if(pamFlag == FALSE){
							progress$inc(1/(sum(lengths(talSites[[1]][]))), detail = paste("Processing TALENs"))
						} else {
							progress$inc(1/(length(unlist(pamSites)) + sum(lengths(talSites[[1]][]))), detail = paste("Processing TALENs"))
						}
						
						if(!is.na(talSites[[1]][[o]][[p]])){
							inExon <- FALSE
							exonID <- exonNumAct
							
							#Kick out target sites where the cut site is not within an exon
							for(r in 1:nrow(exonDF)){
								if((talSites[[1]][[o]][[p]] >= exonDF$exonStart[r]) && (talSites[[1]][[o]][[p]] <= exonDF$exonEnd[r])){
									inExon <- TRUE
									exonID <- exonNumAct
								}
							}
							
							if(inExon){
								#All sites are on the positive strand
								strandId <- "forward"
								
								#Boolean value determining if there is enough sequence context (40 bp upstream, 40 bp downstream) of the cut site
								#in order to run deletion calculation
								contextCondition <- ((talSites[[1]][[o]][[p]] > 40) && 
																		 	(talSites[[1]][[o]][[p]] + 40) <= nchar(geneSeq))
								
								#Check if there's enough context...
								if(!is.na(talSites[[1]][[o]][[p]]) && !is.null(talSites[[1]][[o]][[p]])){
									
									if(contextCondition){
										#Get the sequence context surrounding the cut site
										context <- window(DNAString(geneSeq), talSites[[1]][[o]][[p]], 80)
										
										#Calculate competition for the sequence generated by window()
										slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
										
										#If the score is greater than the threshold, report it
										if(abs(slopeFrame$slopeMH3Plus) >= threshold){
											
											#Get the TALEN target sequence required to target this site
											if(strandId == "forward"){
												talen <- as.character(talSites[[2]][[o]][[1]][[p]])
											} else if(strandId == "complement"){
												talen <- as.character(reverseComplement(DNAString(talSites[[2]][[o]][[1]][[p]])))
											}
											
											formFrame  <- data.frame(Target_Sequence = talen, 
																							 MENTHU_Score = round(abs(slopeFrame$slopeMH3Plus), digits = 2), 
																							 Frame_Shift = slopeFrame$frameShift,
																							 Tool_Type = "TALEN", 
																							 Strand = strandId, 
																							 Exon_ID = exonID, 
																							 Cut_Location = as.integer(talSites[[1]][[o]][[p]], digits = 0))
											
											menthuFrame <- rbind(menthuFrame, formFrame)
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	#colnames(menthuFrame) <- c("Target Sequence", "Score", "Tool", "Strand", "Exon", "Location")
	return(menthuFrame)
}


#' convertToNumeric
#'
#' @param characterStringOfNumbers 
#'
#' @return
#' @export
#'
#' @examples

convertToNumeric <- function(characterStringOfNumbers){
	#Split the string by ','
	splittedString <- unlist(strsplit(characterStringOfNumbers, "[\\, |\\,| ]+"))
	
	#If there is more than one object in the list
	if(length(splittedString) > 1){
		#Find instances where a range was specified
		rangeLocs <- grep("-", splittedString, fixed = TRUE)
		
		#If more than one range
		if(length(rangeLocs) > 1){
			rangeEnds <- strsplit(splittedString[rangeLocs], "-")
			numberList <- splittedString[-rangeLocs]
			
			#Generate a sequence in the range
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


#' Title
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
	if(pamList != ""){
		pamDistList  <- rep(-3, length(pamList))
		custDistList <- strsplit(custList, "[\\, |\\,| ]+")
		distList <- unlist(c(pamDistList, custDistList))
	} else {
		distList <- unlist(strsplit(custList, "[\\, |\\,| ]+"))
	}
	return(as.numeric(distList))
}


#' Title
#'
#' @param pamList 
#' @param custPamList 
#'
#' @return
#' @export
#'
#' @examples
pamStitch <- function(pamList, custPamList){
	if(pamList != ""){
		custList <- strsplit(custPamList, "[\\, |\\,| ]+")
		customList <- unlist(c(pamList, custList))
	} else {
		customList <- unlist(strsplit(custPamList, "[\\, |\\,| ]+"))
	}
	
	return(customList)
}

#' Title
#'
#' @param casList 
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

calculateMENTHUGeneSeqGenBank <- function(casList, cutDistList, wiggle = TRUE, wigRoom = 39, talenList, gbFlag, gbhFlag, 
																					genbankInfo, threshold, firstExon, exonTargetType, exonStuff, progress){
	#Rename the input casList to pamList
	pamList <- casList
	
	#Update progress bar
	progress$inc(0.01, detail = "Scanning for target sites...")
	
	#Get exon sequences and information through getExon
	exon <- getExon(genbankInfo, wiggle = TRUE, wigRoom = 39, gbFlag, exonTargetType, firstExon, exonStuff)
	
	#Get exon indices
	exonInfo <- exon[[1]]
	#Get the exon sequences
	exonSeq  <- exon[[2]]
	#Get the gene sequence
	geneSeq  <- exon[[3]]
	
	#Scan for target sites
	#If the user is using Cas:
	if(length(pamList) > 0){
		#If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
		pamSites <- pamScan(pamList, 
												cutDistList, 
												DNAStringSet(exonSeq), 
												exonList   = exonInfo$exonNum, 
												exonStarts = exonInfo$start, 
												findCut    = TRUE, 
												type       = "cas9", 
												wiggle     = TRUE, 
												wigRoom    = 39)

		#Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
	} else {
		#If the user is NOT using Cas, set pamFlag to FALSE
		pamSites <- 0
		pamFlag <- FALSE
	}
	
	print(talenList)
	#Set a flag to be true if there are TALEN inputs
	talFlag <- unique(talenList != "")
	
	print(talFlag)
	#If there are TALEN inputs
	if(talFlag){
		#Set the range flag to true
		rFlag <- TRUE
		
		#If there are exon inputs
		#Set all exon starts to the exon starts in the input frame
		#Submit talen info to talPal
		talSites <- talPal(DNAStringSet(exonSeq),
											 findcut = TRUE,
											 wiggle = TRUE,
											 wigRoom = 39,
											 range      = rFlag, 
											 armin      = talenList[[1]], 
											 armax      = talenList[[2]], 
											 spamin     = talenList[[3]], 
											 spamax     = talenList[[4]], 
											 exonList   = exonInfo$exonNum,
											 exonStarts = exonInfo$start)
		
		
	} else {
		#If TALENs are not used, set talSites list to empty
		talSites <- ""
	}

	#Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence = as.character(), 
														MENTHU_Score = as.numeric(), 
														Frame_Shift = as.character(), 
														Tool_Type = as.character(), 
														Strand = as.character(), 
														Exon_ID = as.numeric(), 
														Cut_Location = as.integer())
	
	#PAM LEVEL
	if(pamFlag){ #check Cas9 is in use
		
		for(i in 1:length(pamSites[])){
			toolTypeI <- names(pamSites)[i]

			#STRAND LEVEL
			for(j in 1:length(pamSites[[i]][])){
				
				#EXON LEVEL
				if(length(pamSites[[i]][[j]][]) > 0){
					
					for(k in 1:length(pamSites[[i]][[j]][])){
						#Get the number of the exon
						exonNum <- names(pamSites[[i]][[j]])[k]
						exonNumAct <- as.numeric(gsub("Exon ", "", exonNum))
						
						#Only do the calculation if the exon is in the interest list
						
						#SITE LEVEL
						if(length(pamSites[[i]][[j]][[k]][]) > 0){
							
							for(l in 1:length(pamSites[[i]][[j]][[k]][])){
								
								#Make the process bar increment for every potential Cas site
								if(talFlag == FALSE){
									progress$inc(1/(length(unlist(pamSites))), detail = paste("Processing ", names(pamSites[i])))
								} else {
									progress$inc(1/(length(unlist(pamSites)) + sum(lengths(talSites[[1]][]))), detail = paste("Processing ", names(pamSites[i])))
								}

								if(!is.na(pamSites[[i]][[j]][[k]][[l]])){
									inExon <- FALSE
									exonID <- 0
									
									#Kick out target sites where the cut site is not within an exon
									for(m in 1:nrow(exonInfo)){
												
										if((pamSites[[i]][[j]][[k]][[l]] >= exonInfo$start[m]) && 
											 (pamSites[[i]][[j]][[k]][[l]] <= exonInfo$end[m])){
											inExon <- TRUE
											exonID <- exonInfo$exonNum[m]
										}
									}

									if(exonID == exonNumAct){
										#If on the positive strand...
										if(j == 1){
											strandId <- "forward"
											
										} else if (j == 2){
											#If on the negative strand...
											strandId <- "complement" 

										}
										
										#Boolean value determining if there is enough sequence context (40 bp upstream, 40 bp downstream) of the cut site
										#in order to run deletion calculation
										contextCondition <- ((pamSites[[i]][[j]][[k]][[l]] > 40) && 
																				 	(pamSites[[i]][[j]][[k]][[l]] + 40) <= nchar(geneSeq))
										
										#Check if there's enough context...
										if(!is.na(pamSites[[i]][[j]][[k]][[l]]) && !is.null(pamSites[[i]][[j]][[k]][[l]])){
											
											if(contextCondition){
												#Get the sequence context surrounding the cut site
												context <- window(DNAString(geneSeq), pamSites[[i]][[j]][[k]][[l]], 80)
												
												#Calculate competition for the sequence generated by window()
												slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
												
												#If the score is greater than the threshold, report it
												if(abs(slopeFrame$slopeMH3Plus) >= threshold){
													
													#Get the CRISPR target sequence required to target this site
													if(strandId == "forward"){
														crispr <- substr(slopeFrame$seq, (40 + -cutDistList[i] - 19), (40 + -cutDistList[i] + nchar(toolTypeI)))
													} else if(strandId == "complement"){
														crispr <- as.character(substr(Biostrings::reverseComplement(DNAString(slopeFrame$seq)), (40 + -cutDistList[i] - 19),  (40 + -cutDistList[i] + nchar(toolTypeI))))
													}
													
													formFrame  <- data.frame(Target_Sequence = crispr, 
																									 MENTHU_Score = round(abs(slopeFrame$slopeMH3Plus), digits = 2), 
																									 Frame_Shift = slopeFrame$frameShift,
																									 Tool_Type = toolTypeI, 
																									 Strand = strandId, 
																									 Exon_ID = exonID, 
																									 Cut_Location = as.integer(pamSites[[i]][[j]][[k]][[l]], digits = 0))
													menthuFrame <- rbind(menthuFrame, formFrame)
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	if(talFlag){
		#For TALENs
		for(o in 1:length(talSites[[1]][])){
			
			if(length(talSites[[1]][[o]] > 0)){
				#Get the number of the exon
				exonNum <- names(talSites[[2]][o])
				exonNumAct <- as.numeric(gsub("Exon ", "", exonNum))
				#Only do the calculation if the exon is in the interest list
				
				for(p in 1:length(talSites[[1]][[o]])){
					
					#EXON LEVEL
					if(length(talSites[[1]][[o]][[p]][]) > 0){
						if(pamFlag == FALSE){
							progress$inc(1/(sum(lengths(talSites[[1]][]))), detail = paste("Processing TALENs"))
						} else {
							progress$inc(1/(length(unlist(pamSites)) + sum(lengths(talSites[[1]][]))), detail = paste("Processing TALENs"))
						}	
						
						if(!is.na(talSites[[1]][[o]][[p]])){
							inExon <- FALSE
							
							#Kick out target sites where the cut site is not within an exon
							for(r in 1:nrow(exonInfo)){
								if((talSites[[1]][[o]][[p]] >= exonInfo$start[r]) && (talSites[[1]][[o]][[p]] <= exonInfo$end[r])){
									inExon <- TRUE
									exonID <- exonInfo$exonNum[r]
								}
							}
							
							if(inExon){
								if(exonID == exonNumAct){
									
									strandId <- "forward"
									
									#Boolean value determining if there is enough sequence context (40 bp upstream, 40 bp downstream) of the cut site
									#in order to run deletion calculation
									contextCondition <- ((talSites[[1]][[o]][[p]] > 40) && 
																			 	(talSites[[1]][[o]][[p]] + 40) <= nchar(geneSeq))
									
									#Check if there's enough context...
									if(!is.na(talSites[[1]][[o]][[p]]) && !is.null(talSites[[1]][[o]][[p]])){
										
										if(contextCondition){
											#Get the sequence context surrounding the cut site
											context <- window(DNAString(geneSeq), talSites[[1]][[o]][[p]], 80)
											
											#Calculate competition for the sequence generated by window()
											slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
											
											#If the score is greater than the threshold, report it
											if(abs(slopeFrame$slopeMH3Plus) >= threshold){
												
												#Get the CRISPR target sequence required to target this site
												if(strandId == "forward"){
													talen <- as.character(talSites[[2]][[o]][[1]][[p]])
												} else if(strandId == "complement"){
													talen <- as.character(Biostrings::reverseComplement(DNAString(talSites[[2]][[o]][[1]][[p]])))
												}
												
												formFrame  <- data.frame(Target_Sequence = talen, 
																								 MENTHU_Score = round(abs(slopeFrame$slopeMH3Plus), digits = 2), 
																								 Frame_Shift = slopeFrame$frameShift,
																								 Tool_Type = "TALEN", 
																								 Strand = strandId, 
																								 Exon_ID = exonID, 
																								 Cut_Location = as.integer(talSites[[1]][[o]][[p]], digits = 0))
												menthuFrame <- rbind(menthuFrame, formFrame)
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	return(menthuFrame)
}

####Exon Handler for Custom Exon Input####
exonHandler <- function(exonRHandsonTable){
	exonTable <- hot_to_r(exonRHandsonTable)
	exonDF <- exonTable[apply(exonTable, MARGIN = 1, function(x) any(x > 0)),]
	
	return(exonDF)
}



#' window
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
#' This function takes a DNA or RNA sequence as input (along with a parameter specifying the type of sequence) and outputs the complement of the input sequence. E.g., "ATTG" will return "TAAC" if type = "DNA" and "UAAC" if type = "RNA"
#'
#' @param seq A DNA or RNA sequence from which to generate a complement string
#' @param type Default is "DNA"; a DNA sequence can only contain "A", "C", "G", or "T" for the purposes of complement(). The other option is "RNA"; an RNA sequence can only contain "A", "C", "G", or "U" for the purposes of complement().
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
	
	fromVal <- c("A", "C", "G", "T", "a", "c", "g", "t")
	toVal   <- c("T", "G", "C", "A", "t", "g", "c", "a")
	
	
	compSeq <- plyr::mapvalues(unlist(strsplit(compSeq, split = "")),
														 from = fromVal,
														 to   = toVal,
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
#' @param type Default is "DNA"; allowed characters are "A", "C", "G", and "T" (case insensitive). Other option is "RNA"; allowed characters are "A", "C", "G", and "U" (case insensitive.)
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







