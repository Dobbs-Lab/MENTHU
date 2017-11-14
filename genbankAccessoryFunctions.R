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
getExon <- function(genbankInfo, gbFlag, exonTargetType, firstExon, exonStuff) {
	
	require(genbankr)
	
	if(gbFlag){
		# Transform accession string to readable format and parse GenBank data
		#gba <- GBAccession(accession) 
		#gb <- readGenBank(gba)
		seq <- Biostrings:::getSeq(genbankInfo)[[1]] # fetch gene sequence
		exonInfo <- data.frame(slot(exons(genbankInfo),"ranges")) #fetch exon information
	} else {
		seq <- Biostrings:::DNAString(genbankInfo$ORIGIN)
		exonInfo <- getExonLocus(genbankInfo)[, 1:3]
	}
	print(head(exonInfo))
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
			exStart <- exonInfo$start[numExons[i]]
			exEnd <- exonInfo$end[numExons[i]]
			set <- c(set, seq[exStart:exEnd])
			exonSeq <- DNAStringSet(set)
		}
		
		print("here")
	return(list(exonInfo[1:numExons,], exonSeq, seq))
}

#For handling GenBank files that throw errors
wonkyGenBankHandler <- function(gba){
	gbFile <- rentrez:::entrez_fetch(db = "nucleotide", gba, rettype = "gb") 
	geneIn <- unlist(strsplit(gbFile, "\\\n", perl = TRUE))
	
	geneInfo <- formatApe(geneIn)
	return(geneInfo)
}

