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
getExon <- function(genbankInfo, wiggle = TRUE, wigRoom = 39, gbFlag, exonTargetType, firstExon, exonStuff) {
	
	require(genbankr)
	require(Biostrings)
	
	#If the genbank file read properly all the way through...
	if(gbFlag){
		# Transform accession string to readable format and parse GenBank data
		#gba <- GBAccession(accession) 
		#gb <- readGenBank(gba)
		genSeq <- Biostrings:::getSeq(genbankInfo)[[1]] # fetch gene sequence
		exonInfo <- data.frame(slot(exons(genbankInfo),"ranges")) #fetch exon information
	} else {
		#If the genbank file was wonky....
		genSeq <- Biostrings:::DNAString(genbankInfo$ORIGIN) #Get the gene sequence
		exonInfo <- getExonLocus(genbankInfo)[, 1:3] #Get the exon locations
	}
	
	print(head(exonInfo)) #For debugging purposes
	
	#If the user wants to target within a certain percentage of the beginning of the sequence
	if(exonTargetType == 1){
		exonCutoff <- ceiling(nrow(exonInfo) * exonStuff)
		#If the user does not want to count the first exon
		if(firstExon == 0){
			exonMin <- 2
		}	else {
			exonMin <- 1
		}
		exonList <- seq(from = exonMin, to = exonCutoff)
		
		#If the user wants to target within a certain percentage of the end of the sequence
	} else if(exonTargetType == 2){
		exonCutoff <- length(exonInfo)
		exonMin <- exonCutoff - ceiling(exonCutoff * exonStuff)
		exonList <- seq(from = exonMin, to = exonCutoff)
		
		#If the user wants to target a specific list of exons
	} else if(exonTargetType == 3){
		exonList <- convertToNumeric(exonStuff)
		
		#If the user wants to use all exons...
	} else {
		if(firstExon == 0){
			exonList <- seq(from = 2, to = nrow(exonInfo))
		} else {
			exonList <- seq(from = 1, to = nrow(exonInfo))
		}
	}
	
	#print(head(exonList))
		# Variable initialization
		set <- NULL
		exonSeq <- NULL
		
		# Calculation of number of exons to extract from percent parameter
		#numExons <- floor(percent/100*length(exonInfo$start))
		numExons <- exonList
		# Generation of DNAStringSet of exon DNA sequences
		print(length(numExons))
		
		for (i in 1:length(numExons)) {
			print(i)
			#If wiggle = true, include sequence context upstream and downstream of exon so that cut sites whose context runs out of the exon can still be considered
			if(wiggle){
				#Ensure that there is enough wiggle room to add sequence context at beginning of exon (e.g., if exon 1 starts at base 23, there is not 39 bases of wiggle room to add)
				if(exonInfo$start[numExons[i]] - wigRoom < 1){
					exStart <- 1
				} else {
					exStart <- exonInfo$start[numExons[i]] - wigRoom
				}
				
				#Ensure there is enough wiggle room to add sequence context at end of exon (e.g., if exon 10 ends at base 455, and the sequence ends at base 450, there is not 39 bases of wiggle room)
				if(exonInfor$end[numExons[i]] + wigRoom > nchar(seq)){
					exEnd <- nchar(genSeq)
				} else {
					exEnd <- exonInfo$end[numExons[i]] + wigRoom
				}
				set <- c(set, genSeq[exStart:exEnd])
				exonSeq <- DNAStringSet(set)
				
			} else {
				#If no wiggle room...
				exStart <- exonInfo$start[numExons[i]]
				exEnd <- exonInfo$end[numExons[i]]
				set <- c(set, genSeq[exStart:exEnd])
				exonSeq <- DNAStringSet(set)
			}
			
		}
		
		print("here")
	return(list(exonInfo[numExons,], exonSeq, genSeq))
}

#For handling GenBank files that throw errors
wonkyGenBankHandler <- function(gba){
	gbFile <- rentrez:::entrez_fetch(db = "nucleotide", gba, rettype = "gb") 
	geneIn <- unlist(strsplit(gbFile, "\\\n", perl = TRUE))
	
	geneInfo <- formatApe(geneIn)
	return(geneInfo)
}

