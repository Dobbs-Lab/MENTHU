#' getGenbankFile
#'
#' @param accession 
#' @param file 
#' @param deleteTempFile 
#'
#' @return
#' @export
#'
#' @examples

getGenbankFile <- function(accession, file = 'temp.gb', deleteTempFile = TRUE){
	require(curl)
	# Get database from accession format
	#db <- .getDatabaseFromAccession(accession)
	db <- 'nuccore'
	
	# Construct the URL
	baseURL <- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
	rettype <- 'gbwithparts'
	retmode <- 'text'
	tURL    <- paste0(baseURL, "?db=", db, "&id=", accession, "&retmode=", retmode, "&rettype=", rettype)
	
	# Get file from URL using curl
	curl::curl_download(url = tURL, file)
	
	# Read in file contents
	gbContents <- readLines(file)
	
	# Delete temporary file
	if(deleteTempFile){
		file.remove(file)
	}
	
	# Format the file contents
	gbContents <- formatApe(gbContents)
	
	return(gbContents)
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

getExon <- function(genbankInfo, wiggle = TRUE, wiggleRoom = 39, gbFlag, exonTargetType, firstExon, exonStuff) {
	require(Biostrings)
	gbFlag <- FALSE
	
	# If the genbank file read properly all the way through...
	if(gbFlag){
	} else {
		# This is now the handler for all GenBank accessions - genbankr had some issues
		genSeq   <- Biostrings:::DNAString(genbankInfo$ORIGIN) # Get the gene sequence
		exonInfo <- getExonLocus(genbankInfo)
		exonInfo$exonNum <- seq(from = 1, to = nrow(exonInfo))
	}
	
	# Make sure these are numeric to avoid crashes
	exonInfo$start  <- as.numeric(exonInfo$start)
	exonInfo$end    <- as.numeric(exonInfo$end)
	exonInfo$number <- as.numeric(exonInfo$number)
	
	# If there is only one exon in the file, use that exon
	if(nrow(exonInfo) == 1 ){
		exonList <- c(1)
		
	} else {
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
			exonMin    <- exonCutoff - ceiling(exonCutoff * exonStuff)
			exonList   <- seq(from = exonMin, to = exonCutoff)
			
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
	}
	
	# Variable initialization
	set     <- NULL
	
	# Calculation of number of exons to extract from percent parameter
	#numExons <- floor(percent/100*length(exonInfo$start))
	numExons <- exonList
	# Generation of DNAStringSet of exon DNA sequences
	
	for (i in 1:length(numExons)) {
		
		# If wiggle = true, include sequence context upstream and downstream of exon 
		# so that cut sites whose context runs out of the exon can still be considered
		if(wiggle){
			# Ensure that there is enough wiggle room to add sequence context at beginning of exon 
			# (e.g., if exon 1 starts at base 23, there is not 39 bases of wiggle room to add)
			
			if(exonInfo$start[numExons[i]] - wiggleRoom < 1){
				exStart <- 1
			} else {
				exStart <- exonInfo$start[numExons[i]] - wiggleRoom
			}
			
			# Ensure there is enough wiggle room to add sequence context at end of exon 
			# (e.g., if exon 10 ends at base 455, and the sequence ends at base 450, there is not 39 bases of wiggle room)
			if(exonInfo$end[numExons[i]] + wiggleRoom > nchar(as.character(genSeq))){
				exEnd <- nchar(as.character(genSeq))
			} else {
				exEnd <- exonInfo$end[numExons[i]] + wiggleRoom
			}
			
			set <- c(set, as.character(genSeq[exStart:exEnd]))
			
		} else {
			# If no wiggle room...
			exStart <- exonInfo$start[numExons[i]]
			exEnd   <- exonInfo$end[  numExons[i]]
			set     <- c(set, genSeq[exStart:exEnd])
		}
	}
	
	return(list(exonInfo[numExons,], set, as.character(genSeq)))
}


#' wonkyGenBankHandler
#'
#' @param gba 
#'
#' @return
#' @export
#'
#' @examples

wonkyGenBankHandler <- function(gba, apikey = ""){
	require(rentrez)
	# In future releases allowing for high-throughput GenBank retrieval, the UI will allow users to specify an API key
	# to allow for faster returns from NCBI. 
	apikey <- "" # Enable for running locally
	
	# If the user does not supply their own API key:
	if(apikey == ""){
		# Use rentrez API to get genbank flat file coresponding to gba accession
		gbFile <- rentrez:::entrez_fetch(db = "nucleotide", 
																		 gba, 
																		 rettype = "gb", 
																		 retmode = "text")
		# Read the flat file line by line
		geneIn   <- unlist(strsplit(gbFile, "\\\n", perl = TRUE))
		#Pass the readLines from the flat file to formatApe to create a plasmid object
		geneInfo <- formatApe(geneIn) 
		
	} else { # If the user does supply an API key:
		# Use rentrez API to get genbank flat file coresponding to gba accession
		gbFile <- rentrez:::entrez_fetch(db = "nucleotide", 
																		 gba, 
																		 rettype = "gb", 
																		 retmode = "text", 
																		 api_key = apikey)
		# Read the flat file line by line
		geneIn   <- unlist(strsplit(gbFile, "\\\n", perl = TRUE))
		# Pass the readLines from the flat file to formatApe to create a plasmid object
		geneInfo <- formatApe(geneIn)
	}
	
	return(geneInfo)
}

