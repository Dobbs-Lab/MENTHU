#' Title
#'
#' @param id 
#'
#' @return
#' @export
#'
#' @examples
#' 

handleEnsemblInput <- function(id, wiggle = FALSE, wiggleRoom = 39){
	# Get the type (transcript, exon, gene, protein)
	ensType <- getEnsemblIdType(id, check = TRUE)
	
	if(ensType == 'transcript'){
		exonTable   <- processEnsTranscript(id, wiggle, wiggleRoom)
		
	} else if(ensType == 'exon'){
		# Get info about exon, including sequence
		if(wiggle){
			exonTable <- getEnsemblExonSequence(id, wiggleRoom, wiggleRoom)
		} else {
			exonTable <- getEnsemblExonSequence(id)
		}
		
	} else if(ensType == 'protein'){
		# Identify the transcript the protein came from
		parentT <- lookupEnsemblInfo(id)$Parent
		
		# Do transcript stuff to it
		exonTable <- processEnsTranscript(parentT, wiggle, wiggleRoom)
		
	} else if(ensType == 'gene'){
		# Currently, cry and do nothing, but this is alas not a long-term solution
		
	} else if(ensType == 1){
		return(c(1, paste0("Error: The input accession doesn't appear to be an Ensembl accession. ", 
											 "Please check your Enseml input.")))
		
	} else if(ensType == 'unknown'){
		return(c(1, paste0("Error: The type of Ensembl accession can't be determined. ", 
											 "Please check your Ensembl input.")))
		
	} else {
		return(c(1, paste0("Error: This Ensembl accession corresponds to a ", 
											 ensType, 
											 ". Please input an Ensembl transcript, exon, protein, or gene.")))
	}
}

#' getEnsemblGeneFile
#'
#' @param accession 
#' @param file 
#' @param deleteTempFile 
#'
#' @return
#' @export
#'
#' @examples

getEnsemblGeneFile <- function(accession){
	require(httr)
	
	# Construct the URL
	server <- 'http://rest.ensembl.org/sequence/id/'
	
	ext <- paste0(getBaseAccession(accession), "?type=genomic")
	
	# Get file from URL
	r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("text/x-fasta"))
	
	httr::stop_for_status(r)
	
	# Read in file contents
	ensContents <- httr::content(r)
	
	ensContents <- processEnsSequence(ensContents)
	
	return(ensContents)
}


#' getEnsemblExonInfo
#'
#' @param accession 
#'
#' @return
#' @export
#'
#' @examples
#' 

getEnsemblExonInfo <- function(accession){
	require(httr)
	
	# Get the exonic information
	server <- "http://rest.ensembl.org"
	ext    <- paste0("/overlap/id/", getBaseAccession(accession), "?feature=exon")
	
	r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
	
	#httr::stop_for_status(r)
	r <- httr::content(r)
	
	# FOR DEALING WITH ENSEMBL BULLSH*T WHERE THE FIELDS FOR THE EXONS ARE NOT RETURNED IN THE SAME ORDER
	nameKey <- unique(unlist(lapply(r, names)))
	exonInfo <- as.data.frame(setNames(do.call(mapply, c(FUN=c, lapply(r, '[', nameKey))), nameKey), stringsAsFactors = FALSE)
	
	# Identify the "parent" transcript, and filter out any returned exons overlapping this area that are NOT from that transcript
	versionLessId    <- gsub("[.][0-9]+", "", accession, perl = TRUE)
	relevantExonInfo <- exonInfo[which(exonInfo$Parent == versionLessId), ]
	
	return(relevantExonInfo)
}

#' getEnsemblExonInfoG
#'
#' @param accession 
#'
#' @return
#' @export
#'
#' @examples
#' 

getEnsemblExonInfoG <- function(accession){
	require(httr)
	
	# Get the exonic information
	server <- "http://rest.ensembl.org"
	ext    <- paste0("/overlap/id/", getBaseAccession(accession), "?object_type=gene;feature=exon")
	
	s <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
	
	s <- .dealWithEnsemblsBullSh1t(s)
	
	# Get the exonic information
	server <- "http://rest.ensembl.org"
	ext    <- paste0("/overlap/id/", getBaseAccession(accession), "?object_type=gene;feature=transcript")
	
	r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
	
	#httr::stop_for_status(r)
	r <- httr::content(r)
	
	# FOR DEALING WITH ENSEMBL BULLSH*T WHERE THE FIELDS FOR THE EXONS ARE NOT RETURNED IN THE SAME ORDER
	nameKey <- unique(unlist(lapply(r, names)))
	exonInfo <- as.data.frame(setNames(do.call(mapply, c(FUN=c, lapply(r, '[', nameKey))), nameKey), stringsAsFactors = FALSE)
	
	# Identify the "parent" transcript, and filter out any returned exons overlapping this area that are NOT from that transcript
	versionLessId    <- gsub("[.][0-9]+", "", accession, perl = TRUE)
	relevantExonInfo <- exonInfo[which(exonInfo$Parent == versionLessId), ]
	
	return(relevantExonInfo)
}

#' getEnsemblIdType
#'
#' @param id 
#' @param check 
#'
#' @return
#' @export
#'
#' @examples
#' 

getEnsemblIdType   <- function(id, check = FALSE){
	ensemblType <- TRUE
	
	if(check){
		# Check for Ensembl type ID
		ensemblType <- ensemblIdSpecies(id, bool = TRUE)
	}
	
	if(ensemblType){
		# Get the species prefix
		species <- ensemblIdSpecies(id)
		
		# Remove the context
		suffix <- gsub("[0-9.]+", "", gsub(species$Id, "", id))
		
		# Because life is generally unfair, Ensembl and FlyBase (which is where  Ensembl's drosophila 
		# stuff comes from) have different naming schemes. WHICH, BY THE WAY, ARE NOT 
		# EXPLAINED ON ENSEMBL'S STABLE ID PAGE. So THANKS, guys. </sarcasm>
		if(species$Id == "FB"){ # Deal with these differently because it's too much to ask for everyone to
			# have the same naming conventions
			if(       suffix == "TR" | suffix == "tr" | suffix == "tR" | suffix == "Tr"){
				type <- 'transcript'
			} else if(suffix == "GN" | suffix == "gn" | suffix == "gN" | suffix == "Gn"){
				type <- 'gene'
			} else if(suffix == "PP" | suffix == "pp" | suffix == "Pp" | suffix == "Pp"){
				type <- 'protein'
			} #else if(suffix == ){
			
			#} 
			else {
				type <- 'unknown'
			}
			
		} else {
			# For all the not-drosophila entries out there, that conform to these standards:
			if(       suffix == "G"  | suffix == "g"){
				type <- 'gene'
			} else if(suffix == "E"  | suffix == "e"){
				type <- 'exon'
			} else if(suffix == "FM" | suffix == "fm" | suffix == "fM" | suffix == "Fm"){
				type <- 'protein family'
			} else if(suffix == "GT" | suffix == "gt" | suffix == "gT" | suffix == "Gt"){
				type <- 'gene tree'
			} else if(suffix == "P"  | suffix == "p"){
				type <- 'protein'
			} else if(suffix == "R"  | suffix == "r"){
				type <- 'regulatory feature'
			} else if(suffix == "T"  | suffix == "t"){
				type <- 'transcript'
			} else {
				type <- 'unknown'
			}
		}
		return(type)
		
	} else {
		return(1)
	}
}

#' lookupEnsemblInfo
#'
#' @param accession 
#'
#' @return
#' @export
#'
#' @examples
#' 

lookupEnsemblInfo <- function(accession){
	library(httr)
	library(jsonlite)
	library(xml2)
	
	server <- "http://rest.ensembl.org"
	ext <- "/lookup/id/"
	
	r <- httr::GET(paste(server, ext, getBaseAccession(accession), sep = ""), httr::content_type("application/json"))
	
	#httr::stop_for_status(r)
	
	resTable <- data.frame(t(sapply(httr::content(r),c)), stringsAsFactors = FALSE)
	
	return(resTable)
}

#' detectAccessionDB
#' 
#' This function attempts to determine whether an accession comes from Ensembl, GenBank, or RefSeq based 
#' on the accession format and returns 1 if the database cannot be identified based on the information in id. 
#' Note that this only checks the format, and does not check if, for example, the Ensembl prefix actually
#' corresponds to a species.
#'
#' @param id An accession id
#'
#' @return
#' @export
#'
#' @examples
#' 

detectAccessionDB <- function(id){
	
	if(isGenBankFormat(id)){
		return("genbank")
		
	} else if(isRefSeqFormat(id)){
		return("refseq")
		
	} else if(ensemblIdSpecies(id, bool = TRUE)){
		return("ensembl")
		
	} else {
		return("unknown")
		
	}
}

isGenBankFormat <- function(id){
	
	if(grepl("\\b[A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{3}[0-9]{5}|[A-Z]{4}[0-9]{8,10}|[A-Z]{5}[0-9]{7}", id, ignore.case = TRUE, perl = TRUE)){
		return(TRUE)
		
	} else {
		return(FALSE)
		
	}
}

isRefSeqFormat <- function(id){
	if(grepl("[ANWXY][CGTWZMRP][_]", substr(id, 1, 3), ignore.case = TRUE, perl = TRUE)){
		return(TRUE)
		
	} else {
		return(FALSE)
		
	}
}

#' ensemblIdSpecies
#'
#' @param id 
#' @param ensIdList 
#' @param bool 
#'
#' @return
#' @export
#'
#' @examples
#' 

ensemblIdSpecies <- function(id, ensIdList = ensIds, bool = FALSE){
	# Pull out the mouse matches
	if(substr(id, 1, 3) == "MGP"){
		prefix <- unlist(stringr::str_extract_all(id, "[Mm][Gg][Pp][_][A-Za-z0-9]+[_]"))
	} else if(substr(id, 1, 2) == "FB"){
		# Pull out the drosophila matches, which for some inexplicable reason have their own type formatting
		prefix <- "FB"
		
	} else {
		# Get the ID number and the preceding two letters
		suffix <- unlist(stringr::str_extract_all(id, "[A-Za-z][A-Za-z][0-9.]+"))
		suffixTwoLetter <- substr(suffix, 1, 2)
		suffixOneLetter <- substr(suffix, 2, 2)
		
		if(grepl("FM|GT|fm|gt|Fm|Gt|fM|gT", suffixTwoLetter, perl = TRUE)){
			prefix <- gsub(suffix, "", id)
			
		} else if(grepl("E|e|G|g|P|p|R|r|T|t", suffixOneLetter, perl = TRUE)){
			prefix <- gsub(unlist(stringr::str_extract_all(id, "[A-Za-z][0-9.]+")), "", id)
			
		} else {
			
			if(bool){
				return(FALSE)
				
			} else {
				return("")
				
			}
		}
	}
	
	# Match Ensembl id to list of Ensembl Ids
	match <- prefix %in% ensIdList$Id
	
	if(match){
		species <- ensIdList[which(ensIdList$Id == prefix), ]
		
		if(bool){
			return(TRUE)
			
		} else {
			return(species)
			
		}
	} else {
		if(bool){
			return(FALSE)
			
		} else {
			return("unknown")
		}
	}
}



#' isEnsemblUp
#'
#' Function to ping Ensembl's servers and determine if Ensembl is, in fact, up
#' 
#' @return
#' @export
#'
#' @examples
#' 

isEnsemblUp <- function(){
	require(httr)
	require(jsonlite)
	require(xml2)
	
	# This code comes from Ensembl's instructions, pretty much verbatim
	server <- "http://rest.ensembl.org"
	ext    <- "/info/ping?"
	
	r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
	
	httr::stop_for_status(r)
	
	return(sapply(httr::content(r),c))
}


#' getBaseAccession
#' 
#' 

getBaseAccession <- function(accession){
	baseAcc <- strsplit(accession, ".", fixed = TRUE)	
	return(unlist(baseAcc[[1]][1]))
}

#' getEnsemblExonSequence
#'
#' @param accession 
#'
#' @return
#' @export
#'
#' @examples
#' 

getEnsemblExonSequence <- function(accession, strand = 0, exp5 = 0, exp3 = 0){
	require(httr)
	
	# Construct the URL
	server <- 'http://rest.ensembl.org/sequence/id/'
	
	ext <- paste0(getBaseAccession(accession), "?type=genomic;expand_5prime=", exp5, ";expand_3prime=", exp3)
	
	# Get file from URL
	r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("text/x-fasta"))
	
	#httr::stop_for_status(r)
	
	# Read in file contents
	ensContents <- httr::content(r)
	
	ensContents <- processEnsSequence(ensContents)
	ensContents$start   <- as.numeric(ensContents$start)
	ensContents$end     <- as.numeric(ensContents$end)
	
	ensContents$contextStart <- ensContents$start
	ensContents$contextEnd   <- ensContents$end
	ensContents$exp5         <- exp5
	ensContents$exp3         <- exp3
	
	# Fix the start/end info
	if(as.numeric(ensContents$strand) == 1){
		ensContents$start <- ensContents$start + exp5
		ensContents$end   <- ensContents$end   - exp3
		
	} else if(as.numeric(ensContents$strand) == -1){
		ensContents$start <- ensContents$start + exp5
		ensContents$end   <- ensContents$end   - exp3 
	} 
	
	return(ensContents)
}

#' processEnsTranscript
#'
#' @param accession An Ensembl TRANSCRIPT ID (last letter in the prefix is a "T")
#' @param wiggle If true, additional bases will be added to the examined exons so that target sites that would otherwise run off the end of the exon can still be targeted
#' @param wiggleRoom The number of additional bases of sequence to consider if wiggle is true
#'
#' @return
#' @export
#'
#' @examples
#' 

processEnsTranscript <- function(id, wiggle = TRUE, wiggleRoom = 39){
	# Get info about exons in transcript
	exonTable <- getEnsemblExonInfo(id)
	
	# Get the strand
	strand <- exonTable$strand[1]
	
	# Get the exons and sequences from ensembl
	if(wiggle){
		exonSeqs <- lapply(exonTable$exon_id, getEnsemblExonSequence, strand, wiggleRoom, wiggleRoom)
	} else {
		exonSeqs <- lapply(as.character(exonTable$exon_id), getEnsemblExonSequence, strand)
	}
	
	# Put the retrieved exon sequences into a data frame
	exonSeqCols        <- names(exonSeqs[[1]])
	exonSeqs           <- as.data.frame(matrix(unlist(exonSeqs), 
																						 nrow  = length(exonSeqs), 
																						 ncol  = length(unlist(exonSeqs[[1]])), 
																						 byrow = TRUE), 
																			stringsAsFactors = FALSE)
	colnames(exonSeqs) <- exonSeqCols
	
	# merge the exonTable and exonSeqs tables 
	exonTable <- merge(exonTable, exonSeqs, by = c("start", "end"))
	
	
	# Get the exon information
	#exonTable <- getEnsemblExonInfo(accession)
	
	# Get the sequence
	#transSeq  <- getEnsemblGeneFile(accession)
	
	# Figure out where the gene starts and stops
	#geneStart <- as.numeric(transSeq$start)
	#geneEnd   <- as.numeric(transSeq$end)
	
	# Map the genomic coordinates to the gene sequence
	#if(transSeq$strand == "-1"){ # If the gene is on the reverse strand
	#	exonStartT <- abs(as.numeric(unlist(exonTable$end))   - geneEnd   - 1)
	#	exonEndT   <- abs(as.numeric(unlist(exonTable$start)) - geneEnd   - 1)	
	
	#} else if(transSeq$strand == "1"){ # If the gene is on the forward strand
	#	exonStartT <- abs(as.numeric(unlist(exonTable$start)) - geneStart + 1)
	#	exonEndT   <- abs(as.numeric(unlist(exonTable$end))   - geneStart + 1)
	#}	
	
	# Create a data frame with the re-indexed exons
	#newExonTable <- data.frame(
	#	parent    = unlist(exonTable$Parent),
	#	exonId    = unlist(exonTable$exon_id),
	#	exonNum   = as.numeric(unlist(exonTable$rank)),
	#	exonStart = exonStartT,
	#	exonEnd   = exonEndT,
	#	stringsAsFactors = FALSE
	#)
	
	# If there should be wiggle room, add the additional sequence if possible
	#if(wiggle){
	#	exonSeqs <- sapply(1:nrow(newExonTable), function(x) substr(transSeq$sequence[1], 
	#																															(if((newExonTable$exonStart[x] - wiggleRoom) >= 1){
	#																																newExonTable$exonStart[x] - wiggleRoom
	#																															} else {
	#																																1
	#																															}), 
	#																															(if((newExonTable$exonEnd[x]   + wiggleRoom) <= nchar(transSeq$sequence)){
	#																																newExonTable$exonEnd[x]   + wiggleRoom
	#																															} else {
	#																																nchar(transSeq$sequence)
	#																															})))
	
	#} else {
	#	exonSeqs <- sapply(1:nrow(newExonTable), function(x) substr(transSeq$sequence[1], newExonTable$exonStart[x], newExonTable$exonEnd[x]))
	#}
	
	# Add the exon sequences to the exon table
	#newExonTable$exonSeqs <- exonSeqs
	
	return(exonTable)
}


#' processEnsProtein
#'
#' @param accession An Ensembl PROTEIN ID (last letter in the prefix is a "P")
#' @param wiggle If true, additional bases will be added to the examined exons so that target sites that would otherwise run off the end of the exon can still be targeted
#' @param wiggleRoom The number of additional bases of sequence to consider if wiggle is true
#'
#' @return
#' @export
#'
#' @examples
#' 

processEnsProtein <- function(accession, wiggle = TRUE, wigRoom = 39){
	# Get the parent transcript
	infoTable <- lookupEnsemblInfo(accession)
	parentT   <- infoTable$Parent
	
	# Process as a transcript
	proteinInfo <- processEnsTranscript(parentT, wiggle, wigRoom)
	
	return(proteinInfo)
}

processEnsGene <- function(){
	
}


processEnsExon <- function(accession, wiggle = TRUE, wiggleRoom = 39){
	if(wiggle){
		exonSeq <- getExonSequence(accession, wiggleRoom, wiggleRoom)
		
	} else {
		exonSeq <- getExonSequence(accession)
		
	}
	
}

#' Title
#'
#' @param inSeq 
#'
#' @return
#' @export
#'
#' @examples
#' 

processEnsSequence <- function(inSeq){
	# Split the fasta header and sequence body
	header       <- strsplit(inSeq, "\\n")[[1]][1]
	sequenceBody <- gsub("\\n", "", gsub(header, "", inSeq))
	
	# Process the information in the fasta header
	splitList <- unlist(strsplit(header, " |:"))
	id        <- gsub(">", "", splitList[1])
	version   <- splitList[3]
	chromNum  <- splitList[4]
	seqStart  <- splitList[5]
	seqEnd    <- splitList[6]
	strand    <- splitList[7]
	
	# Put the information into a datatable
	seqTable <- data.frame(id = id,
												 genomeVersion = version,
												 chromosome    = chromNum,
												 strand        = strand,
												 start         = seqStart,
												 end           = seqEnd,
												 sequence      = sequenceBody,
												 stringsAsFactors = FALSE)
	
	return(seqTable)
}


#' .dealWithEnsemblsBullSh1t
#' 
#' For some godforsaken reason, when multiple items are returned from Ensembl (e.g., a gene with two or more associated transcripts),
#' the information fields are not always in the same order. WHHHHHHHHYYYYYYYYYYYYYY. This function makes tables by putting that info in the same order.
#'
#' @param ensContents 
#'
#' @return
#' @export
#'
#' @examples
#' 

.dealWithEnsemblsBullSh1t <- function(ensContents){
	
	r <- httr::content(ensContents)
	
	# FOR DEALING WITH ENSEMBL BULLSH*T WHERE THE FIELDS FOR THE EXONS ARE NOT RETURNED IN THE SAME ORDER
	nameKey <- unique(unlist(lapply(r, names)))
	nameKey <- nameKey[which(nameKey != "description")]
	
	# Drop the description field if it is present
	ensTable <- as.data.frame(setNames(do.call(mapply, c(FUN=c, lapply(r, '[', nameKey))), nameKey), stringsAsFactors = FALSE)
	
	return(ensTable)
}

#' getChromosomeInfo
#'
#' @return
#' @export
#'
#' @examples
#' 

getChromosomeInfo <- function(species){
	# Get the chromosome information for homo_sapiens
	server <- "http://rest.ensembl.org"
	ext    <- "/info/assembly/homo_sapiens?"
	
	r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
	
	# Get the info from the retrieval
	struc <- sapply(httr::content(r), c)
	
	# Format Ensembl return 
	strucChrome <- struc$top_level_region[union(union(which(lapply(struc$top_level_region, "[[", 1) == "chromosome"), 
																										which(lapply(struc$top_level_region, "[[", 2) == "chromosome")),  
																							which(lapply(struc$top_level_region, "[[", 3) == "chromosome"))]
	
	# Generate a data frame of the chromosome data
	chromeTable <- as.data.frame(matrix(unlist(strucChrome), ncol = 3))
	
	# Convert chromosome length to numeric
	chromeTable <- as.data.frame(dplyr::bind_rows(strucChrome), stringsAsFactors = FALSE)
	chromeTable$length <- as.numeric(chromeTable$length)
	
	# Order by chromosome length
	chromeTable <- chromeTable[order(-chromeTable$length), ]
	
	return(chromeTable)
}


processChromosome <- function(species, chromosome, chromeTable, wiggle = FALSE, wiggleRoom = 0){
	# Generate server info for next retrieval
	server     <- "http://rest.ensembl.org" 
	ext        <- "/overlap/region/"
	
	# Get the info for the current chromosome
	chromInfo <- chromeTable[which(chromeTable$name == chromosome), ]
	
	# Generate a table to contain information about the regions
	regions <- seq(from = 1, to = chromInfo$length, by = 5000000)
	
	# Fill in the end
	regionTable           <- data.frame(regionStart = regions, stringsAsFactors = FALSE)
	regionTable$regionEnd <- c(seq(from = 5000000, to = chromInfo$length, by = 5000000), chromInfo$length)
	
	returnData <- list()
	
	# Get the exons within the region
	for(i in 1:nrow(regionTable)){
		print(i)
		
		serverCall <- paste0(server, 
												 ext, 
												 species, "/", 
												 chromosome, ":", 
												 format(regionTable$regionStart[i], scientific = FALSE), "-", 
												 format(regionTable$regionEnd[i],   scientific = FALSE), 
												 "?feature=exon;")
		
		r2 <- httr::GET(serverCall, httr::content_type("application/json"))
		r3 <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r2)))
		returnData <- c(returnData, r3)
	}
	
	# Generate unique list of exons
	idList <- unique(returnData$exon_id)
	
	
	if(wiggle){
		exonSeqs <- getEnsemblExonSequence(idList[1], exp5 = wiggleRoom, exp3 = wiggleRoom)
		# Get the exon sequences for all exons
		for(j in 2:length(idList)){
			tExon    <- getEnsemblExonSequence(idList[j], exp5 = wiggleRoom, exp3 = wiggleRoom)
			
			exonSeqs <- rbind(exonSeqs, tExon)
		}
		
	} else {
		exonSeqs <- getEnsemblExonSequence(idList[1], exp5 = wiggleRoom, exp3 = wiggleRoom)
		
		# Get the exon sequences for all exons
		for(j in 2:length(idList)){
			tExon    <- getEnsemblExonSequence(idList[j])
			
			exonSeqs <- rbind(exonSeqs, tExon)
		}
	}
	
	
	return(exonSeqs)
	
}