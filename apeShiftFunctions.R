#apeShift functions
#' formatApe
#' Function to take the input of a GenBank .gb or Plasmid Editor .ape flat file and format into a data frame
#' Note that this function currently treats unknown upper and lower boundaries as known - e.g.,
#' a location specified as '<1..121' indicates that the feature starts before the first sequenced base, and continues
#' to base 121 (inclusive). However, this is treated as the feature consisting of bases 1-121.
#' Also does not currently support join(complement(), complement()) files, but DOES support
#' complement(join()) type files.
#'
#' @param apeContents
#'
#' @return
#' @export
#'
#' @examples
#'

formatApe <- function(apeContents){
	#Identify where the origin sequence starts
	#featureLoc     <- grep("^FEATURES",       apeContents)
	
	#Get location of all top-level fields
	#topFieldList  <- grep("^[^\\s]",          apeContents[1:(featureLoc - 1)], perl = TRUE)
	#subFieldList  <- grep("^[\\s]{1,3}[\\S]", apeContents[1:(featureLoc - 1)], perl = TRUE)
	#fullFieldList <- sort(c(topFieldList, subFieldList))
	
	#fieldTable <- data.frame(fieldName = character(length(fullFieldList)), fieldValue = character(length(fullFieldList)))
	#for(i in 1:length(fullFieldList)){
	#  matchTerm <- strsplit(apeContents[fullFieldList[i]], " ")[[1]][1]
	#  if(matchTerm != ""){
	#    fieldTable$fieldName[i] <- matchTerm
	#  }
	#  matchTerm <- stripWhiteSpace(matchTerm)
	#}
	
	#Get fields for the class, and assign their values
	locus      <- sub("\\s+",     "", gsub("LOCUS",      "", grep("^LOCUS",        apeContents, value = TRUE)))
	#name      <- strsplit(locus, "\\s{2,}", perl = TRUE)[[1]][2]
	definition <- sub("\\s+",     "", gsub("DEFINITION", "", grep("^DEFINITION",   apeContents, value = TRUE)))
	accession  <- sub("\\s+",     "", gsub("ACCESSION",  "", grep("^ACCESSION",    apeContents, value = TRUE)))
	version    <- sub("\\s+",     "", gsub("VERSION",    "", grep("^VERSION",      apeContents, value = TRUE)))
	source     <- sub("\\s+",     "", gsub("SOURCE",     "", grep("^SOURCE",       apeContents, value = TRUE)))
	organism   <- sub("\\s+",     "", gsub("ORGANISM",   "", grep("^\\s*ORGANISM", apeContents, value = TRUE)))
	#reference  <- sub("\\s+",     "", gsub("REFERENCE",  "", grep("^REFERENCE",    apeContents, value = TRUE)))
	#authors    <- sub("\\s+",     "", gsub("AUTHORS",    "", grep("^\\s*AUTHORS",  apeContents, value = TRUE)))
	#title      <- sub("\\s+",     "", gsub("TITLE",      "", grep("^\\s*TITLE",    apeContents, value = TRUE)))
	
	#Get the reference(s)
	#references <- grep()
	
	
	comments   <- gsub("\\s{2,}", "", gsub("COMMENT",    "", grep("^COMMENT",      apeContents, value = TRUE)))
	comments   <- comments[which(comments != "")]
	
	
	#featureIndexStart <- grep("^FEATURES", apeContents, value = FALSE)
	featureIndexStop  <- grep("^ORIGIN",   apeContents, value = FALSE)
	
	#Get the indices of each feature
	featureIndices  <- grep("^\\s+[a-zA-Z0-9_\\-\\']+\\s+[a-zA-Z\\>\\<]*\\({0,1}[0-9]+[\\.][\\.][\\>\\<]{0,1}[0-9]+[\\>\\<]{0,1}\\){0,1}", apeContents)
	featureIndices1 <- grep("^\\s+[a-zA-Z0-9_\\-\\']+\\s+complement", apeContents)
	featureIndices2 <- grep("^\\s+[a-zA-Z0-9_\\-\\']+\\s+join", apeContents)
	featureIndices  <- unique(sort(c(featureIndices, featureIndices1, featureIndices2)))
	
	#Get the sequence
	seqLines <- grep("[0-9]+", apeContents[featureIndexStop:length(apeContents)], value = TRUE, perl = TRUE)
	seqStart <- as.numeric(strsplit(gsub("\\s{2,}", "", seqLines[1]), " ")[[1]][1])
	seq <- gsub("\\s+",   "", seqLines)
	seq <- gsub("[0-9]+", "", seq)
	seq <- toupper(paste(seq, collapse = ''))
	
	featList <- list()
	
	#Get all values associated with each feature
	for(i in 1:length(featureIndices)){
		curFeat <- strsplit(gsub("^\\s{2,}", "", apeContents[featureIndices[i]]), "\\s{2,}", perl = TRUE)
		attType <- curFeat[[1]][1]
		
		#Get the index/indices of the feature
		if(i < length(featureIndices)){
			indexStop <- featureIndices[i + 1] - 1
		} else {
			indexStop <- featureIndexStop - 1
		}
		
		featureString <- paste(apeContents[featureIndices[i]:(indexStop - 1)], collapse = "")
		indexString   <- strsplit(featureString, "\\/", fixed = FALSE)[[1]][1]
		indexString   <- gsub("^\\s*[0-9]+", "", indexString, perl = TRUE)
		indexString   <- gsub("[a-zA-Z]*", "", indexString, perl = TRUE)
		indexString   <- gsub("\\'", "", indexString, perl = TRUE)
		indexString   <- gsub("\\s*", "", indexString, perl = TRUE)
		indexString   <- gsub("\\(", "", indexString,  perl = TRUE)
		indexString   <- gsub("\\)", "", indexString,  perl = TRUE)
		indexString   <- gsub("\\>", "", indexString,  perl = TRUE)
		indexString   <- gsub("\\<", "", indexString,  perl = TRUE)
		
		#If the indices of the feature are complicated by a join
		if(grepl("join\\(", curFeat, ignore.case = TRUE)){
			#If on the complement strand
			if(grepl("complement\\(", curFeat, ignore.case = TRUE)){
				orientation <- "complement"
				#Get the string containing the join info
				#joinString <- grep("complement\\([a-zA-Z0-9\\,\\.\\s\\>\\<]+\\)", perl = TRUE, ignore.case = TRUE, value = TRUE)
			} else {
				#String is in default orientation
				orientation <- "default"
				#Get the string containing the join info
				#joinString <- grep("join\\([0-9\\,\\.\\s\\>\\<]+\\)", joinS, perl = TRUE, ignore.case = TRUE, value = TRUE)
			}
			#Do a bunch of cleaning steps to get the start and end joining indices into a list
			#Remove "complement("
			#jS <- gsub("complement\\(", "", joinString)
			#Remove "join("
			#jS <- gsub("join", "", jS)
			#Remove any lone "("
			#jS <- gsub("\\(", "", jS)
			#Remove ending ")"
			#jS <- gsub("\\)", "", jS)
			#Get rid of newline characters
			#jS <- gsub("\\\n", "", jS)
			#Remove any ">" or "<" characters
			#jS <- gsub("\\>", "", jS, perl = TRUE)
			#jS <- gsub("\\<", "", jS, perl = TRUE)
			#Split into start and stop indices
			jS <- strsplit(indexString, split = ",", fixed = TRUE)
			jSStart <-  as.numeric(sapply(strsplit(unlist(jS), split = "\\.\\.", fixed = FALSE), "[", 1))
			jSEnd   <-  as.numeric(sapply(strsplit(unlist(jS), split = "\\.\\.", fixed = FALSE), "[", 2))
			
			seqFeat <- paste(substr(rep(seq, length(jSStart)), jSStart, jSEnd), collapse = "")
			
			#Input NAs for startIn and stopIn
			startIn <- NA
			stopIn  <- NA
			
		} else {
			if(grepl("complement", curFeat, ignore.case = TRUE)){
				orientation <- "complement"
				startIn <- as.numeric(                strsplit(gsub("complement\\(", "", gsub("\\<", "", gsub("\\>", "", indexString, perl = TRUE), perl = TRUE)), "\\.\\.", perl = TRUE)[[1]][1])
				stopIn  <- as.numeric(gsub("\\)", "", strsplit(gsub("complement\\(", "", gsub("\\<", "", gsub("\\>", "", indexString, perl = TRUE), perl = TRUE)), "\\.\\.", perl = TRUE)[[1]][2]))
				
			} else {
				orientation <- "default"
				startIn <- as.numeric(strsplit(gsub("\\<", "", gsub("\\>", "", indexString, perl = TRUE)), "\\.\\.", perl = TRUE)[[1]][1])
				stopIn  <- as.numeric(strsplit(gsub("\\<", "", gsub("\\>", "", indexString, perl = TRUE)), "\\.\\.", perl = TRUE)[[1]][2])
			}
			jSStart <- NA
			jSEnd   <- NA
			
			#Get the corresponding sequence
			seqFeat <- substr(seq, startIn, stopIn)
		}
		
		
		if(i < length(featureIndices)){
			featValues <- getFeatureValues(apeContents[featureIndices[i]:(featureIndices[i + 1] - 1)])
		} else {
			featValues <- getFeatureValues(apeContents[featureIndices[i]:(featureIndexStop - 1)])
		}
		
		featValues[1, 1] <- "feature_type"
		featValues$value[1] <- gsub("\"", "", attType)
		featValues[2, 1] <- "featStart"
		featValues$value[2] <- startIn
		featValues[3, 1] <- "featEnd"
		featValues$value[3] <- stopIn
		featValues[4, 1] <- "orientation"
		featValues$value[4] <- orientation
		featValues[5, 1] <- "joinStart"
		featValues$value[5] <- list(jSStart)
		featValues[6, 1] <- "joinStop"
		featValues$value[6] <- list(jSEnd)
		featValues[7, 1] <- "genomicContext"
		featValues$value[7] <- toupper(as.character(seqFeat))
		featValues[8, 1] <- "featureSequence"
		featValues$value[8] <- (if(orientation == "complement"){toupper(reverseComplement(as.character(seqFeat)))} else {toupper(as.character(seqFeat))})
		featList[[i]] <- featValues
	}
	
	
	
	ape <- list(locus, definition, accession, version, source, organism, comments, featList, seqStart, seq)
	names(ape) <- list("LOCUS", "DEFINITION", "ACCESSION", "VERSION", "SOURCE", "ORGANISM", "COMMENT", "FEATURES", "seqStart", "ORIGIN")
	
	class(ape) <- "apePlasmid"
	
	return(ape)
}



getFeatures <- function(plasmid, qual = NULL, qualValue = NULL){
	if(!is.null(qual) || !is.null(qualValue)){
		plasRet <- list()
		
		if(!is.null(qual)) {
			plasRet <- plasmid$FEATURES[grepl(paste(qual, collapse = "|"), lapply(plasmid$FEATURES, "[[", 1))]
		}
		
		if(!is.null(qualValue)) {
			if(length(plasRet) > 0){
				
				plasQualVal <- plasmid$FEATURES[grepl(paste(qualValue, collapse = "|"), lapply(plasmid$FEATURES, "[[", 2))]
				plasRet <- append(plasRet, plasQualVal)
				
			} else {
				plasRet <- plasmid$FEATURES[grepl(paste(qualValue, collapse = "|"), lapply(plasmid$FEATURES, "[[", 2))]
				
			}
		}
		return(plasRet)
		
	} else {
		
		return(plasmid$FEATURES)
	}
}

#' getFeatureValues
#'
#' @param featureLines
#'
#' @return
#' @export
#'
#' @examples
#'

getFeatureValues <- function(featureLines){
	require(stringr)
	#Get all the lines that start a qualifier
	valueLines <- grep("/", featureLines, perl = TRUE)
	#Create a data frame to hold the qualifiers
	df <- data.frame(qualifier = character(length(valueLines) + 8), value = c(length(valueLines) + 8), stringsAsFactors = FALSE)
	#Write dynamic qualifiers after the mandatory first five qualifiers (feature_type, featStart, featEnd, orientation, join)
	row <- 9
	
	df$value <- NA
	
	#For each qualifier
	for(i in valueLines){
		#Get the name of the qualifier
		valName <- gsub("\\s+/", "", strsplit(featureLines[i], "=", fixed = TRUE)[[1]][1])
		
		#Determine if the qualifier spans multiple lines, and stitch the lines together if so
		if(stringr:::str_count(featureLines[i], "\"") == 1){
			searchLine <- paste0(featureLines[i], featureLines[i + 1])
		} else {
			searchLine <- featureLines[i]
		}
		
		#Get the value for the qualifier
		valVal <- gsub("\\s+", " ", strsplit(searchLine, "=", fixed = TRUE)[[1]][2])
		valVal <- gsub("\"", "", valVal)
		
		#Put in data frame
		df[row, 1] <- valName
		df[row, 2] <- valVal
		#Advance the row counter
		row <- row + 1
	}
	
	return(df)
}



getExonLocus <- function(gene){
	geneF <- getFeatures(gene)
	exonListI <- which(sapply(sapply(geneF, "[", 2), "[", 1) == "exon")
	geneExons <- geneF[exonListI]
	exonTable <- data.frame(start = numeric(length(exonListI)), 
													end   = numeric(length(exonListI)),
													width = numeric(length(exonListI)),
													type  = character(length(exonListI)),
													orientation = character(length(exonListI)),
													number = character(length(exonListI)),
													stringsAsFactors = FALSE)
	
	for(i in 1:length(exonListI)){
		exonTable[i, 1] <- geneExons[[i]][2, 2]
		exonTable[i, 2] <- geneExons[[i]][3, 2]
		exonTable[i, 3] <- as.numeric(geneExons[[i]][3, 2]) - as.numeric(geneExons[[i]][2, 2])
		exonTable[i, 4] <- geneExons[[i]][1, 2]
		exonTable[i, 5] <- geneExons[[i]][4, 2]
		exonTable[i, 6] <- geneExons[[i]][which(geneExons[[i]]$qualifier == "number"),2]
	}
	
	return(exonTable)
}




