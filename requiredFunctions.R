####Required Function#######

calculateMENTHUGeneSeq <- function(casList, wiggle = TRUE, wigRoom = 39, geneSeq, threshold, exonDF, progress, armin, armax, spamin, spamax){
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
			pamSites <- pamScan(pamList, exonSeqs, exonStarts = exonDF$exonStart, findCut = TRUE, type = "cas9", wiggle = wiggle, wigRoom = wigRoom)
		} else {
			pamSites <- pamScan(pamList, exonSeqs, exonStarts = NULL, findCut = TRUE, type = "cas9", wiggle = wiggle, wigRoom = wigRoom)
		}
		
		#Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
	} else {
		
		#If the user is NOT using Cas, set pamFlag to FALSE
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
												 exonStarts = exonDF$exonStart)
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
												 exonStarts = NULL)
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
												crispr <- substr(slopeFrame$seq, 24, (43 + nchar(toolTypeI)))
											} else if(strandId == "complement"){
												crispr <- as.character(substr(reverseComplement(DNAString(slopeFrame$seq)), 24, 43 + nchar(toolTypeI)))
											}
											
											#Create data frame of the current results
											formFrame  <- data.frame(Target_Sequence = crispr, 
																							 MENTHU_Score = round(abs(slopeFrame$slopeMH3Plus), digits = 2), 
																							 Frame_Shift = slopeFrame$frameShift,
																							 Tool_Type = toolTypeI, 
																							 Strand = strandId, 
																							 Exon_ID = exonNum, 
																							 Cut_Location = as.integer(pamSites[[i]][[j]][[k]][[l]], digits = 0))
											
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
					print(paste0("P = ", p))
					print(paste0("talSites NOP = ", talSites[[1]][[o]][[p]]))
					
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
								print(exonDF$exonStart[r])
								print(talSites[[1]][[o]][[p]])
								if((talSites[[1]][[o]][[p]] >= exonDF$exonStart[r]) && (talSites[[1]][[o]][[p]] <= exonDF$exonEnd[r])){
									inExon <- TRUE
									exonID <- exonNumAct
								}
							}
							
							print(inExon)
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

convertToNumeric <- function(characterStringOfNumbers){
	splittedString <- unlist(strsplit(characterStringOfNumbers, ","))
	if(length(splittedString) > 1){
		rangeLocs <- grep("-", splittedString, fixed = TRUE)
		if(length(rangeLocs) > 1){
			rangeEnds <- strsplit(splittedString[rangeLocs], "-")
			numberList <- splittedString[-rangeLocs]
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

calculateMENTHUGeneSeqGenBank <- function(casList, talenList, gbFlag, gbhFlag, genbankInfo, threshold, firstExon, exonTargetType, exonStuff, progress){
	#Rename the input casList to pamList
	pamList <- casList
	
	#Subset geneSeq to exon sequences
	#exonSeqs <- substring(rep(as.character(geneSeq), nrow(exonDF)), exonDF$exonStart, exonDF$exonEnd)
	
	#Update progress bar
	progress$inc(0.1, detail = "Scanning for target sites...")
	
	#if(gbFlag){
	#	geneSeq <- as.character(info@sequence)
	#} else {
	#	geneSeq <- as.character(info$ORIGIN)
	#}

	exon     <- getExon(genbankInfo, gbFlag, exonTargetType, firstExon, exonStuff)
	#Get exon indices
	exonInfo <- exon[[1]]
	#Get the exon sequences
	exonSeq  <- exon[[2]]
	#Get the gene sequence
	geneSeq  <- exon[[3]]
	
	#Generate the list of exons that should be checked
	
	print("here 1")
	#Scan for target sites
	pamSites <- pamScan(pamList, DNAStringSet(exonSeq), findCut = TRUE, type = "cas9")	
	#pamSites <- pamScan(pamList, DNAStringSet(exonSeqs), findCut = TRUE, type = "cas9")
	
	#if(length(talenList > 0)){
	#targetList, findcut = FALSE, range = FALSE, armin = 15, armax = 18, spamin = 14, spamax = 16, exonStarts = NULL
	#	talenSites <- talPal(targetList, 
	#											 findcut = TRUE, 
	#											 range   = TRUE, 
	#											 armin   = talenList[1], 
	#											 armax   = talenList[2], 
	#											 spamin  = talenList[3], 
	#											 spamax  = talenList[4],
	#											 exonStarts = NULL)
	#}
	
	#Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence = as.character(), 
														MENTHU_Score = as.numeric(), 
														Frame_Shift = as.character(), 
														Tool_Type = as.character(), 
														Strand = as.character(), 
														Exon_ID = as.numeric(), 
														Cut_Location = as.integer())
	print(pamSites)
	#PAM LEVEL
	for(i in 1:length(pamSites[])){
		print(i)
		toolTypeI <- names(pamSites)[i]
		progress$inc(1/(length(pamSites[])), detail = paste("Processing ", names(pamSites[i])))
		#print(paste0("i = ", i)) #For testing purposes
		#print(nchar(geneSeq))    #For testing purposes
		#STRAND LEVEL
		for(j in 1:length(pamSites[[i]][])){
			print(j)
			#EXON LEVEL
			if(length(pamSites[[i]][[j]][]) > 0){
				
				
				for(k in 1:length(pamSites[[i]][[j]][])){
					print(k)
					
					#Get the number of the exon
					exonNum <- names(pamSites[[i]][[j]])[k]
					exonNumAct <- as.numeric(gsub("Exon ", "", exonNum))
					#Only do the calculation if the exon is in the interest list
					#if(exonNumAct %in% exonList){
					#SITE LEVEL
					if(length(pamSites[[i]][[j]][[k]][]) > 0){
						
						for(l in 1:length(pamSites[[i]][[j]][[k]][])){
							print(l)
							if(!is.na(pamSites[[i]][[j]][[k]][[l]])){
								#print(l) #For testing purposes
								#inExon <- TRUE
								exonID <- exonNumAct
								#Kick out target sites where the cut site is not within an exon
								for(m in 1:nrow(exonInfo)){
									if((pamSites[[i]][[j]][[k]][[l]] >= exonInfo$start[m]) && (pamSites[[i]][[j]][[k]][[l]] <= exonInfo$end[m])){
										inExon <- TRUE
										exonID <- exonNumAct
									}
								}
								
								#if(inExon){
								#If on the positive strand...
								if(j == 1){
									strandId <- "forward"
									#print(strandId) #For testing purposes
									#print(pamSites[[i]][[j]][[k]][[l]]) #For testing purposes
									
								} else if (j == 2){
									#If on the negative strand...
									strandId <- "complement" 
									#print(strandId) #For testing purposes
									#print(paste0(pamSites[[i]][[j]][[k]][[l]])) #For testing purposes
								}
								
								#Boolean value determining if there is enough sequence context (40 bp upstream, 40 bp downstream) of the cut site
								#in order to run deletion calculation
								contextCondition <- ((pamSites[[i]][[j]][[k]][[l]] > 40) && 
																		 	(pamSites[[i]][[j]][[k]][[l]] + 40) <= nchar(geneSeq))
								
								#Check if there's enough context...
								if(!is.na(pamSites[[i]][[j]][[k]][[l]]) && !is.null(pamSites[[i]][[j]][[k]][[l]])){
									#print(contextCondition) #For testing purposes
									
									if(contextCondition){
										#Get the sequence context surrounding the cut site
										context <- window(DNAString(geneSeq), pamSites[[i]][[j]][[k]][[l]], 80)
										#print(context) #For testing purposes
										
										#Calculate competition for the sequence generated by window()
										slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
										#print(abs(slopeFrame$slopeMH3Plus)) #For testing purposes
										
										#If the score is greater than the threshold, report it
										if(abs(slopeFrame$slopeMH3Plus) >= threshold){
											
											#Get the CRISPR target sequence required to target this site
											if(strandId == "forward"){
												crispr <- substr(slopeFrame$seq, 24, (43 + nchar(toolTypeI)))
											} else if(strandId == "complement"){
												crispr <- as.character(substr(reverseComplement(DNAString(slopeFrame$seq)), 24, 43 + nchar(toolTypeI)))
											}
											
											#print(crispr) #For testing purposes
											#print(round(abs(slopeFrame$slopeMH3Plus), digits = 2)) #For testing purposes
											
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
	
	
	#For TALENs
	for(i in 1:length(talSites[])){
		print(i)
		toolTypeI <- names(talSites)[i]
		progress$inc(1/(length(talSites[])), detail = paste("Processing ", names(talSites[i])))
		#print(paste0("i = ", i)) #For testing purposes
		#print(nchar(geneSeq))    #For testing purposes
		#STRAND LEVEL
		for(j in 1:length(talSites[[i]][])){
			print(j)
			#EXON LEVEL
			if(length(talSites[[i]][[j]][]) > 0){
				
				
				for(k in 1:length(talSites[[i]][[j]][])){
					print(k)
					
					#Get the number of the exon
					exonNum <- names(talSites[[i]][[j]])[k]
					exonNumAct <- as.numeric(gsub("Exon ", "", exonNum))
					#Only do the calculation if the exon is in the interest list
					#if(exonNumAct %in% exonList){
					#SITE LEVEL
					if(length(talSites[[i]][[j]][[k]][]) > 0){
						
						for(l in 1:length(talSites[[i]][[j]][[k]][])){
							print(l)
							if(!is.na(talSites[[i]][[j]][[k]][[l]])){
								#print(l) #For testing purposes
								#inExon <- TRUE
								exonID <- exonNumAct
								#Kick out target sites where the cut site is not within an exon
								for(m in 1:nrow(exonInfo)){
									if((talSites[[i]][[j]][[k]][[l]] >= exonInfo$start[m]) && (talSites[[i]][[j]][[k]][[l]] <= exonInfo$end[m])){
										inExon <- TRUE
										exonID <- exonNumAct
									}
								}
								
								#if(inExon){
								#If on the positive strand...
								if(j == 1){
									strandId <- "forward"
									#print(strandId) #For testing purposes
									#print(talSites[[i]][[j]][[k]][[l]]) #For testing purposes
									
								} else if (j == 2){
									#If on the negative strand...
									strandId <- "complement" 
									#print(strandId) #For testing purposes
									#print(paste0(talSites[[i]][[j]][[k]][[l]])) #For testing purposes
								}
								
								#Boolean value determining if there is enough sequence context (40 bp upstream, 40 bp downstream) of the cut site
								#in order to run deletion calculation
								contextCondition <- ((talSites[[i]][[j]][[k]][[l]] > 40) && 
																		 	(talSites[[i]][[j]][[k]][[l]] + 40) <= nchar(geneSeq))
								
								#Check if there's enough context...
								if(!is.na(talSites[[i]][[j]][[k]][[l]]) && !is.null(talSites[[i]][[j]][[k]][[l]])){
									#print(contextCondition) #For testing purposes
									
									if(contextCondition){
										#Get the sequence context surrounding the cut site
										context <- window(DNAString(geneSeq), talSites[[i]][[j]][[k]][[l]], 80)
										#print(context) #For testing purposes
										
										#Calculate competition for the sequence generated by window()
										slopeFrame <- calculateSlopeCompetition(as.character(context), cutSite = 40, weight = 20, top = 10)
										#print(abs(slopeFrame$slopeMH3Plus)) #For testing purposes
										
										#If the score is greater than the threshold, report it
										if(abs(slopeFrame$slopeMH3Plus) >= threshold){
											
											#Get the CRISPR target sequence required to target this site
											if(strandId == "forward"){
												crispr <- substr(slopeFrame$seq, 24, (43 + nchar(toolTypeI)))
											} else if(strandId == "complement"){
												crispr <- as.character(substr(reverseComplement(DNAString(slopeFrame$seq)), 24, 43 + nchar(toolTypeI)))
											}
											
											#print(crispr) #For testing purposes
											#print(round(abs(slopeFrame$slopeMH3Plus), digits = 2)) #For testing purposes
											
											formFrame  <- data.frame(Target_Sequence = crispr, 
																							 MENTHU_Score = round(abs(slopeFrame$slopeMH3Plus), digits = 2), 
																							 Frame_Shift = slopeFrame$frameShift,
																							 Tool_Type = toolTypeI, 
																							 Strand = strandId, 
																							 Exon_ID = exonID, 
																							 Cut_Location = as.integer(talSites[[i]][[j]][[k]][[l]], digits = 0))
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
	print(menthuFrame) #For testing purposes
	#colnames(menthuFrame) <- c("Target Sequence", "Score", "Tool", "Strand", "Exon", "Location")
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







