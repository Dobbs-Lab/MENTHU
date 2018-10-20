# Required packages
library(shiny)
library(shinyjs)
library(shinyTable)
library(rhandsontable)
library(shinyIncubator)
library(stringr)
library(stringi)
library(Biostrings)
library(rentrez)
library(rlist)
library(DT)
library(plyr)

# Required supporting files
source("apeShiftFunctions.R")
source("genbankAccessoryFunctions.R")
source("menthu2.0AccessoryFunctions.R")
source("required2.0Functions_1.R")
source("targetAccessoryFunctions2.0.R")
<<<<<<< HEAD
=======
source("ensemblAccessoryFunctions.R")
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"

shinyServer(function(input, output, session){
	########################################################
	###################Global Variables#####################
	########################################################
	
	# Storage for final results so that it is accessible to download handler and submit buttons
	results   <<- 0
	
	# Empty data frame for displaying when users are inputting exon information
	dfEmpty   <<- data.frame(Exon_Num         = rep(0, 5),
												   exonStart        = rep(0, 5), 
												   exonEnd          = rep(0, 5), 
												   stringsAsFactors = FALSE)
	
	# Example data frame for when the paste gene seq example is selected
	dfExample <<- data.frame(Exon_Num         = c(1,     2,   3,   4,   5),
													 exonStart        = c(1,   101, 201, 301, 401),
													 exonEnd          = c(100, 200, 300, 400, 500),
													 stringsAsFactors = FALSE)
	
	# Reactive values
	rValues <- reactiveValues(downloadF          = FALSE,   # Flag for displaying download button (for copy/paste gene seq)
														downloadFGB        = FALSE,   # Flag for displaying download button (for genbank)
														geneSeqResultsFlag = FALSE,   # Flag for displaying gene seq results table
														genbankResultsFlag = FALSE,   # Flag for displaying genbank results table
														rhFrame            = dfEmpty, # Slot to hold exon information data frame
<<<<<<< HEAD
														resetVal           = FALSE,   #  for if the reset button has been clicked
														geneSeqError       = 0)       # Error messages for geneSeq submission
=======
														resetVal           = FALSE,   # For if the reset button has been clicked
														geneSeqError       = 0,       # Error messages for geneSeq submission
														downloadFEns       = FALSE,   # Flag for displaying download button (for Ensembl)
														ensemblResultsFlag = FALSE,   # Flag for displaying Ensembl results table
														validExonListFlag  = TRUE
	)
	
	# Load Ensembl ID table
	ensIds <<- readRDS("2018-09-21_ensIds.RDS")
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
	
	# Flag for displaying example table when clicking example geneSeq link
	geneSeqExampleFlag <<- FALSE
	
	# Data file containing gene database; for use with pre-computed genes tab
	#geneDF  <<-  readRDS(file = "testGeneDB.rds")
	
	########################################################
	##################Validation Checks#####################
	########################################################
	
	####Make sure GenBank/RefSeq ID is properly formatted####
	validRefSeqGBGenbankId <- reactive({
		if(input$genbankId != ""){
			#Let RefSeq RNA accessions through
			if(stringr::str_detect(input$genbankId, regex("^(NM|NR|XM|XR)_[0-9]{6}", ignore_case = TRUE))){
	
				#Catch RefSeq protein accesssions
			} else if(stringr::str_detect(input$genbankId, regex("^(AP|NP|YP|XP|WP)_[0-9]{6}", ignore_case = TRUE))){
				validate(
					need(1 == 2,
							 "Error: This ID matches the RefSeq protein accession format. Please submit an accession corresponding to a DNA sequence."))
				
				# Catch ginormous genomic region files; disabled for running locally
			} else if(stringr::str_detect(input$genbankId, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}", ignore_case = TRUE))){
				#validate(
				#	need(1 == 2,
				#			 paste0("Error: This ID matches a RefSeq complete genomic molecule, incomplete genomic region, contig, ", 
				#			 			 "scaffold, or complete genome. We do not currently support any of these reference types due to ", 
				#			 			 "issues surrounding exon identification. You can use a GenBank or RefSeq nucleotide entry ", 
				#			 			 "corresponding to your gene of interest."))
				#)
			
				# Catch GenBank protein accessions
			} else if(stringr::str_detect(input$genbankId, regex("^[A-Z]{3}[0-9]{5}", ignore_case = TRUE))){
				validate(
					need(1 == 2,
							 "Error: This ID matches GenBank protein accession format. Please submit an accession corresponding to a DNA sequence."))
				
				# Catch anything else not conforming to input types
			} else {
				validate(
					need((((stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{2}[0-9]{6}",                  ignore_case = TRUE)))  | 
								 (stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{1}[0-9]{5}",                  ignore_case = TRUE)))) |
								 (stringr::str_detect(input$genbankId, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}\\.[0-9]", ignore_case = TRUE)))) |
							 	 (stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{4}[0-9]{8,10}",               ignore_case = TRUE))),
							 paste0("Error: The entered ID does not match a known Genbank NUCLEOTIDE or RefSeq NUCLEOTIDE ID format.", 
							 			 "Please check your submitted accession to make sure you are using a NUCLEOTIDE entry.")))
			}
		} 
	})
	
<<<<<<< HEAD
=======
	####Validate Ensembl accession format####
	validEnsemblId <- reactive({
		# Check that the input is not empty
		if((input$inputType == 3) & input$ensemblId != ""){
			# If gene input...
			if(getEnsemblIdType(input$ensemblId) == "gene"){
				shiny::validate(
					need(1 == 2,
							 paste0("Error: The input Ensembl accession appears to be an Ensembl gene. ", 
							 			 "Ensembl genes can have several associated transcripts; please use the Ensembl transcript ID instead of the gene.")
					)
				)
				
				# If NOT gene input...	
			} else {
				shiny::validate(
					need(ensemblIdSpecies(input$ensemblId, ensIds, bool = TRUE),
							 "Error: Input Ensembl ID does not match a known Ensembl species or Ensembl ID format. Please check your input ID."
							 
					)
				)
			}
		}
	})
	
	# Check if the Ensembl ID actually exists
	ensemblIdExists <- reactive({
		if((input$inputType == 3) & input$ensemblId != ""){
			if(is.null(validEnsemblId())){
				
				# If we can make contact with Ensembl
				if(isEnsemblUp()){
					shiny::validate(
						need(!grepl("not found", lookupEnsemblInfo(input$ensemblId)),
								 "Error: The input Ensembl ID was not found in Ensembl's database. Please check your input ID."
						)
					)
					
				# If we can't contact Ensembl
				} else {
					shiny::validate(
						need(1 == 2,
								 "Warning: Ensembl is not responding to our requests. Please try again in a few minutes."
						)
					)
				}
			}
		}
	})
	
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
	####Validate copy/paste sequence input####
	validGeneSeq <- reactive({
		#If the input type is copy/paste gene seq and text has been entered
		if((input$geneSeq   != "") && 
			 (input$inputType == 2)){
			validate(
				#Check for fasta format
				need(!stringr::str_detect(input$geneSeq, "[>]"), 
						 "Error: Input DNA sequence appears to be in FASTA format. Please paste your sequence without the fasta header."),
				
				#Check for not DNA input
				need(!stringr::str_detect(input$geneSeq, "[^ACGTacgt0-9\\s\\n]"), 
						 paste0("Error: Input DNA sequence contains non-standard nucleotides. ", 
						 			 "Allowed nucleotides are A, C, G, and T.")),
				
				# Prevent users from blowing up the server; disabled  for running locally
				#need(nchar(input$geneSeq) < 5000, 
				#		 paste0("The DNA sequence has >5000 nucleotides. For sequences of this size,",
				#		 			 " please use the local version of MENTHU, which can be accessed via the 'Tools and Downloads' tab.")),
				
				# Prevent users from submitting too short sequences
				need(nchar(input$geneSeq) > 80,
						 paste0("The DNA sequence has <80 nucleotides. MENTHU requires at least 40 nucleotides upstream", 
						 			 "and 40 nucleotides downstream of the DSB site in order to properly calculate the MENTHU score."))#,
			
				#Prevent users from blowing up the server
				#need(nchar(input$geneSeq) < 5000, 
				#		 paste0("The DNA sequence has >5000 nucleotides. For sequences of this size, ",
				#		 			  "please use the local version of MENTHU, which can be accessed via the 'Tools and Downloads' tab."))
			)
			
		} else if(input$inputType == 2) {
			#prevent running on empty submission
			validate(
				need(input$geneSeq != "", "")
			)
		}
	})
	
	# Valid threshold - make sure the threshold is not negative
	#validThreshold <- reactive({
	#	if(!is.null(input$threshold)){
	#		validate(
	#			need(input$threshold >= 0,
	#					 "Error: Threshold value must be non-negative.")
	#		)
	#	}	
	#})
	
	# Valid PAM - make sure user selects at least one target type
	validPAM <- reactive({
		validate(
			#need((length(input$casType) > 0) | (input$talenOp == 1) | input$customCutOpt == 1,
			need((length(input$casType) > 0) | input$customCutOpt == 1,
					 #"Error: No nuclease selected. Please select a Cas type and/or a custom PAM scheme and/or the TALEN input option.")
					 "Error: No nuclease selected. Please select a Cas type or choose to use a custom PAM scheme.")
		)
	})
	
	#Valid TALEN inputs
	#validTalen <- reactive({
	#	validate(
			#need(input$armin   <= input$armax,
			#		 "Error: Maximum TALEN arm length must be greater than or equal to minimum TALEN arm length."),
			
			#need(input$spamin  <= input$spamax,
			#		 "Error: Maximum spacer length must be greater than or equal to minimum spacer length."),
			
			#need((input$spamin >= 14),
			#		 "Error: Minimum spacer length is 14 nucleotides."),
			
			#need((input$spamax <= 16),
			#		 "Error: Maximum spacer length is 16 nucleotides."),
			
	#		need((input$armin  >= 15),
	#				 "Error: Minimum TALEN arm length is 15 nucleotides."),
			
	#		need((input$armax  <= 18),
	#				 "Error: Maximum TALEN arm length is 18 nucleotides.")
			
	#	)
	#})
	
	# Valid custom PAM inputs
	validCustomPam <- reactive({
		# If the custom PAM is selected, make sure it only has IUPAC nucleotides or separating characters
		if(input$customPamSeq != ""){
			if(stringr::str_detect(input$customPamSeq, "[^ACGTRYSWKMBDHVNacgtryswkmbdhvn\\s,]")){
				validate(
					need(1 == 2,
							 paste0("Error: Non-allowed characters detected. Only standard nucleotide (A, C, G, T), IUPAC extended ", 
							 			 "nucleotide symbols (R, Y, S, W, K, M, B, D, H, V, N), and separating characters (spaces and commas) allowed."))
				)
			} else {
				# Check for all 'N' PAMs, single nucleotide PAMs, or PAMs consisting solely of Ns and a single nucleotide; disabled for running locally
				#disallowed      <- "^[N]+[ACGTRYSWKMBDHV]{1}$|^[N]+[ACGTRYSWKMBDHV]{1}[N]+$|^[ACGTRYSWKMBDHVN]{1}$|^[N]+$"
				#disallowedCheck <- grepl(disallowed, pamStitch("", input$customPamSeq), ignore.case = TRUE, perl = TRUE)
				#disallowedExist <- is.element(TRUE, unlist(disallowedCheck))
				
				#if(disallowedExist){
				#	validate(
				#		need(1 == 2,
				#				 paste0("Error: Due to computational limitations, we currently do not accept custom PAMs consisting", 
				#				 			 "solely of 'N', single nucleotide PAMs, or PAMs consisting of solely of 'N's and a single nucleotide."))
				#	)
				#}
			}
			
		}
	})
	
	# Valid custom cut sites
	validCustomCutSites <- reactive({
		if(input$cutSite != ""){
			validate(
				need(!stringr::str_detect(input$cutSite, "[^0-9\\s,-]"),
						 "Error: Only numbers, dashes, commas, and spaces allowed.")
			)
		}
	})
	
	# Check that the number of custom PAMs matches the number of custom cutSites
	validMatchCustomInputLength <- reactive({
		if((input$customPamSeq != "") && 
			 (input$customCutOpt != "")){
			validate(
				need(length(as.character(distStitch("", input$cutSite))) == length(pamStitch("", input$customPamSeq)),
				paste0("Error: The number of custom PAM sequences does not match the number of custom DSB sites. ", 
							 "Please make sure that each PAM has one distance to its cut site specified."))
			)
		}
	})
	
	validExonList <- reactive({
		if(input$exonTargetType == 3 && input$exonTargetList != ""){
			if(!rValues$validExonListFlag){
				shiny::validate(
					need(1 == 2,
							 paste0("Error: Some or all of the exons in this list are not found in the input accession. ",
							 			  "Please check the exon numbering of your accession (ESPECIALLY if you're using Ensembl), ",
							 			  "as GenBank, RefSeq, and Ensembl may have inconsistent exon numbering. "))
				)
			}
		}
	})
	
	########################################################
	##############PRINT VALIDATION RESULTS##################
	########################################################
	output$validgenbankid              <- renderText({
		#validGenBankId()
		validRefSeqGBGenbankId()
	})
	
	output$genbankidexists             <- renderText({
		
	})
	
	output$validgeneseq                <- renderText({
		validGeneSeq()
	})
	
	#output$validthreshold              <- renderText({
	#	validThreshold()
	#})
	
	output$validpam                    <- renderText({
		validPAM()
	})
	
	#output$validtalen                  <- renderText({
	#	validTalen()
	#})
	
	output$validcustompam              <- renderText({
		validCustomPam()
	})
	
	output$validcustomcutsites         <- renderText({
		validCustomCutSites()
	})
	
	output$validmatchcustominputlength <- renderText({
		validMatchCustomInputLength()
	})
	
	output$geneseqerrors               <- renderText({
		geneSeqErrors()
	})
	
	output$validexonlist <- renderText({
		validExonList()
	})
	
	#output$validexoninfo <- renderText({
	#	validExonInfo()
	#})
	
	########################################################
	#####################UI Rendering#######################
	########################################################
	
	geneSeqErrors <- reactive({
		rValues$geneSeqError
		if(rValues$geneSeqError == 0){
			""
		} else if(rValues$geneSeqError == 1){
			"No targets corresponding to the PAM(s) detected."
		}
	})
	
	# Render function for input table; DO NOT TOUCH ANY OF THIS - it breaks when you look at it sideways
	renderTable <- reactive({
		exonOutDF <- NULL
		
		if(!is.null(input$exonInfo)){
			exonOutDF <- hot_to_r(input$exonInfo)
			
		} else if(!is.null(isolate(rValues$rhFrame))){
			exonOutDF <- isolate(rValues$rhFrame)
		}
		
		if(!is.null(exonOutDF)){
			
			rValues$rhFrame <- exonOutDF
		}
		
		return(exonOutDF)
	}) %>% debounce(1000)
	
	# Render exon input table; DO NOT TOUCH ANY OF THIS - it took literal months to get this working properly
	output$exonInfo <- renderRHandsontable({
		# Action buttons that this function is dependent on
		input$reset
		input$exampleGeneSeq
		
		if(isolate(rValues$resetVal)){
			exonOutDF <- rValues$rhFrame
		} else {
			exonOutDF <- renderTable()
			rValues$resetVal <- FALSE
		}
		
		if(!is.null(exonOutDF)){
			rhandsontable(exonOutDF) %>%
				hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>% 
				hot_col("Exon_Num",  format = "0") %>%
				hot_col("exonStart", format = "0") %>%
				hot_col("exonEnd",   format = "0")
		}
	})
	

	
	# Download button for copy/paste results
	output$downOutGS <- renderUI({
		if(rValues$downloadF){
			downloadButton("downRes", "Download Results")
		} else {
			""
		}
	})
	
	# Download button for genbank results
	output$downOutGB <- renderUI({
		if(rValues$downloadFGB){
			downloadButton("downRes", "Download Results")
		} else {
			""
		}
	})
	
<<<<<<< HEAD
=======
	# Download button for Ensembl results
	output$downOutEns <- renderUI({
		if(rValues$downloadFEns){
			downloadButton("downRes", "Download Results")
		} else {
			""
		}
	})
	
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
	# Output function for copy/paste results
	output$geneSeqResults <- renderUI({
		if(rValues$geneSeqResultsFlag){
			output$placeholder <- DT::renderDataTable(results, 
																								options  = list(scrollX = TRUE), 
																								rownames = FALSE,
																								escape   = FALSE)
			DT::dataTableOutput("placeholder")
		} else {
			""
		}
		
	})
	
	# Output function for genbank results
	output$genbankResults <- renderUI({
		if(rValues$genbankResultsFlag){
			output$placeholder <- DT::renderDataTable(results, 
																								options  = list(scrollX = TRUE), 
																								rownames = FALSE,
																								escape   = FALSE)
			DT::dataTableOutput("placeholder")
		} else {
		  ""
		}
	})
	
<<<<<<< HEAD
	# Download handler
	output$downRes <- downloadHandler(
		filename = function(){
			#Name file in the form "YYYY-MM-DD_HH-MM-SS_targets.csv
			paste(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_targets.csv")},
		
		content = function(file){
			
			resOut     <- results
			resOut[,1] <- gsub("<strong>|</strong>", "", results[,1], ignore.case = TRUE, perl = TRUE)
			#write.csv( , file, row.names = FALSE, append = FALSE)
			write.table(resOut, file, row.names = FALSE, append = TRUE, quote = FALSE, sep = ",")
		}
	)
	
	# Pre-computed genes functions
	
	# Gene Selection UI Output
	#output$pcGeneSelect <- renderUI({
	#	if(!is.null(input$pcSpeciesSelect)){
	#		if(input$pcSpeciesSelect == "human"){
	#			selectInput("pcGene", 
	#									label = "Please select your gene of interest: ",
	#									choices = c("Not yet supported" = 0))
	#		} else if(input$pcSpeciesSelect == "drerio"){
	#			selectInput("pcGene", 
	#									label = "Please select your gene of interest: ",
	#									choices = c("gene1",
	#															"gene2"))
	#		}
	#	}
	#})
	
	#Clickable Plot Output
	#output$genePlot <- renderPlot({
	#	subPlot <<- geneDF[which(geneDF$Species == input$pcSpeciesSelect & geneDF$GeneId == input$pcGene), ]
	#	barplot(subPlot$Score, names.arg = subPlot$Index, xlab = "Position", ylab = "MENTHU Score", col=ifelse(subPlot$Score >= 1.5,"blue", "gray"))
	#})
	
	#output$siteInfo <- renderText({
	#	paste0("This site is at position ", subPlot$Index[input$plot_single_click$x], " and has a MENTHUv2.0 score of ", subPlot$Score[input$plot_single_click$x])
	#})
	
	#output$plotHoverInfo <- renderPrint({
	#	cat("")
	#})
=======
	# Output function for Ensembl results
	output$ensemblResults <- renderUI({
		if(rValues$ensemblResultsFlag){
			output$placeholder <- DT::renderDataTable(results, 
																								options  = list(scrollX = TRUE),
																								rownames = FALSE,
																								escape   = FALSE)
			DT::dataTableOutput("placeholder")
		} else {
			""
		}
	})
	
	# Download handler
	output$downRes <- downloadHandler(
		filename = function(){
			#Name file in the form "YYYY-MM-DD_HH-MM-SS_InputID_targets.csv
			if(input$inputType == 1){
				# For GenBank/Refseq
				paste(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_", input$genbankId, "_MENTHU_targets.csv")
				
			} else if(input$inputType == 3){
				# For Ensembl
				paste(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_", input$ensemblId, "_MENTHU_targets.csv")
				
			} else {
				# For copy/paste
				paste(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_custom_seq_MENTHU_targets.csv")
			}
		},
		
		content = function(file){
			
			resOut     <- results
			# Remove HTML comments
			resOut[,1] <- gsub("<strong>|</strong>", "", results[,1], ignore.case = TRUE, perl = TRUE)
			# Output the file
			write.table(resOut, file, row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = ",")
		}
	)
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
	
	########################################################
	#################Submission Handling####################
	########################################################
	
	observeEvent(input$genbankSubmit,{
		#talFlag <- 2
		#Run checks for okay PAM/TALEN input
		#if(input$talenOp == 1){
		#	if(is.null(validTalen())){
		#		talFlag <- 1 #Input TALEN is good
		#	} else {
		#		talFlag <- 2 #Problems with input TALEN
		#	}
		#} else {
			talFlag <- 0 # TALENs not selected
		#}
		
		cusPamFlag <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) && is.null(validCustomCutSites())){
				cusPamFlag <- 1
			} else {
				cusPamFlag <- 2
			}
		} else {
			cusPamFlag <- 0
		}
		
		lenMatch <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validMatchCustomInputLength())){
				lenMatch <- 1
			} else {
				lenMatch <- 2
			}
		} else {
			lenMatch <- 0
		}
		
		# Prevent the whole shebang from running without okay inputs
		if(is.null(validRefSeqGBGenbankId()) && # Check Genbank ID is okay
			 #is.null(validThreshold()) &&         # Check threshold is okay
			 is.null(validPAM()) &&               # Check that one of the input options is selected
			 (cusPamFlag != 2) &&                 # Check if custom PAMs are used, and if they are okay
			 (talFlag != 2) &&                    # Check if TALENs are used, and if they are okay
			 (lenMatch != 2)){                    
			
			# Create a new progress object
			progress <- shiny::Progress$new()
			
			# Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			# Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			talenList <- ""
			
			#if(input$talenOp == 1){
			#	if(input$spacer == 0){
			#		spamin <- 14
			#		spamax <- 14
			#	} else if(input$spacer == 1){
			#		spamin <- 16
			#		spamax <- 16
			#	} else {
			#		spamin <- 14
			#		spamax <- 16
			#	}
			#	talenList <- list(input$armin, input$armax, spamin, spamax)
			#}
			
			progress$set(detail = "Retrieving GenBank entry...", value = 0.1)
			
			# Try to pull genbank entry associated with accession
			# Get the GenBank sequence with exon/intron information
			gba <- input$genbankId
			
			endOfTry <<- FALSE
			gbFlag   <<- FALSE
			gbhFlag  <<- FALSE
			
			#Try to retrieve the Genbank file 
			if(endOfTry == FALSE){
				tryCatch({
					info <- suppressWarnings(wonkyGenBankHandler(gba))
					endOfTry <<- TRUE
					gbhFlag  <<- TRUE  # Flag to indicate wonkyGenBankHandler was required 
					gbFlag   <<- FALSE # Flag to indicate readGenBank failed #deprecated; to remove
				}, error = function(err){
				}
				)
			}
			
			# Figure out which exons to target
			if(endOfTry){
				if(input$exonTargetType == 0){
					exonStuff <- 1
					
				} else	if(input$exonTargetType == 1){
					exonStuff <- input$exonBegPercentage / 100
					
				} else if(input$exonTargetType == 2){
					exonStuff <- input$exonEndPercentage / 100
					
				} else if(input$exonTargetType == 3){
					exonStuff <- input$exonTargetList	
				}
				
				# Handle cut distances and PAMs for input to calculateMENTHUGeneSeqGenBank
				if(input$customCutOpt == 1){ # If customs pams in use
					suppressWarnings(if(!is.null(input$casType)){ # If pre-made PAMs in use
						pams         <- pamStitch( input$casType, input$customPamSeq)
						cutDistances <- distStitch(input$casType, input$cutSite)
						
					} else { # If no pre-made PAMs
						print(paste0("Custom PAM Seq: ", input$customPamSeq))
						pams         <- pamStitch( "", input$customPamSeq)
						cutDistances <- distStitch("", input$cutSite)
					}
					)
				} else { # If no custom pams
					pams         <- input$casType
					cutDistances <- rep(-3, length(input$casType))
				}
					
				#if(input$talenOp == 1){
				#	talFlag <- TRUE
				#} else {
					talFlag <- FALSE
				#}
				#Calculate the MENTHU score
				#results <<- calculateMENTHUGeneSeqGenBank(pams, cutDistances, wiggle = TRUE, wigRoom = 39, talenList, gbFlag, gbhFlag, talFlag,
				#																					info, input$threshold, input$firstExon, input$exonTargetType, exonStuff, 
				#																					progress, input$scoreScheme)
				#results <<- calculateMENTHUGeneSeqGenBank(pams, cutDistances, wiggle = TRUE, wigRoom = 39, talenList, gbFlag, gbhFlag, 
				#																					info, input$threshold, input$firstExon, input$exonTargetType, exonStuff, progress)
				stuff    <- calculateMENTHUGeneSeqGenBank(pams, cutDistances, wiggle = TRUE, wigRoom = 39, talenList, gbFlag, gbhFlag, 
																									info, input$firstExon, input$exonTargetType, exonStuff, progress)
				

				# Statistics on the number of targets detected
				output$genbankhits <- renderUI({
					HTML(paste(
						paste0("Number of target sites detected: ",                                             stuff[[2]]),
						paste0("Number of target sites with sufficient sequence context to calculate score: ",  stuff[[3]]),
						paste0("Number of target sites with 3bp microhomology arms within 5bp of each other: ", stuff[[4]]),
						paste0("Number of target sites with score above 1.5 threshold: ",                       stuff[[5]]),
						paste0("Number of target sites satisfying 3bp mh and threshold constraints: ",          stuff[[6]]),
						sep = "<br>"
					))
				})
				
				results <<- stuff[[1]]
				
				# Order the result table from largest menthuScore to smallest, and drop 0s
				results <<- results[which(results$MENTHU_Score > 0), ]
				results <<- results[order(-results$MENTHU_Score), ]

				# Set the flag to display genbank results
				rValues$genbankResultsFlag <- TRUE
				
				# Set the download button flag to true to render download button visible
				rValues$downloadFGB <- TRUE
				
			} else {
				# Output error if no genbank file is found
				output$genbankIdOutcome <- renderText(paste0("Error: Accession '", input$genbankId, "' was not found in database."))
			}
		}
	})
	
<<<<<<< HEAD
=======
	observeEvent(input$ensemblSubmit,{
		timeX <- Sys.time()
		#talFlag <- 2
		#Run checks for okay PAM/TALEN input
		#if(input$talenOp == 1){
		#	if(is.null(validTalen())){
		#		talFlag <- 1 #Input TALEN is good
		#	} else {
		#		talFlag <- 2 #Problems with input TALEN
		#	}
		#} else {
		talFlag <- 0 # TALENs not selected
		#}
		
		cusPamFlag <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) && is.null(validCustomCutSites())){
				cusPamFlag <- 1
			} else {
				cusPamFlag <- 2
			}
		} else {
			cusPamFlag <- 0
		}
		
		lenMatch <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validMatchCustomInputLength())){
				lenMatch <- 1
			} else {
				lenMatch <- 2
			}
		} else {
			lenMatch <- 0
		}
		
		# Prevent the whole shebang from running without okay inputs
		if(is.null(validEnsemblId())  && # Check Ensembl ID is okay
			 is.null(ensemblIdExists()) && # Make sure the ID actually exists
			 #is.null(validThreshold()) && # Check threshold is okay
			 is.null(validPAM())        && # Check that one of the input options is selected
			 (cusPamFlag != 2)          && # Check if custom PAMs are used, and if they are okay
			 (talFlag    != 2)          && # Check if TALENs are used, and if they are okay
			 (lenMatch   != 2)){                    
			
			# Create a new progress object
			progress <- shiny::Progress$new()
			
			# Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			# Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			talenList <- ""
			
			#if(input$talenOp == 1){
			#	if(input$spacer == 0){
			#		spamin <- 14
			#		spamax <- 14
			#	} else if(input$spacer == 1){
			#		spamin <- 16
			#		spamax <- 16
			#	} else {
			#		spamin <- 14
			#		spamax <- 16
			#	}
			#	talenList <- list(input$armin, input$armax, spamin, spamax)
			#}
			
			progress$set(detail = "Retrieving Ensembl entry...", value = 0.1)
			
			# Make sure Ensembl is up and responsive
			if(isEnsemblUp()){
				
				ensemblInfo <- handleEnsemblInput(input$ensemblId, wiggle = TRUE, wiggleRoom = 39)
				
				progress$set(detail = "Processing Ensembl entry...", value = 0.1)
				
				# If the entry is NOT an exon
				if(getEnsemblIdType(input$ensemblId, check = TRUE) != "exon"){
					
					ensemblInfo$rank <- as.numeric(ensemblInfo$rank)
					
					maxR <- max(ensemblInfo$rank)
					
					if(input$firstExon == 1 | max(ensemblInfo$rank == 1)){
						minR <- 1
					} else {
						minR <- 2
					}
					
					# Figure out which exons to target
					if(input$exonTargetType == 0){
						exonStuff <- seq(minR, maxR)
						
					} else if(input$exonTargetType == 1){
						exonStuff <- seq(minR,   ceiling(minR + (maxR * input$exonBegPercentage / 100)))
						
					} else if(input$exonTargetType == 2){
						exonStuff <- seq(max(minR, floor(maxR - (maxR * input$exonEndPercentage / 100))), maxR)
						
					} else if(input$exonTargetType == 3){
						exonStuff <- unique(c(minR, convertToNumeric(input$exonTargetList)))
						
						# Make sure that all of the exons in the list are actually in what is retrieved from Ensembl
						if(all(exonStuff %in% ensemblInfo$rank)){
							rValues$validExonListFlag <- TRUE
							
						} else {
							# Flag that the ensembl exon list is invalid
							rValues$validExonListFlag <- FALSE
							
						}
					}
					
				} else {
					# For when the entry is an exon
					exonStuff <- 1
					
				}
				
				if(rValues$validExonListFlag){
				# Handle cut distances and PAMs for input to calculateMENTHUEnsembl
				if(input$customCutOpt == 1){ # If custom PAMs in use
					suppressWarnings(if(!is.null(input$casType)){ # If pre-made PAMs are also in use
						pams         <- pamStitch( input$casType, input$customPamSeq)
						cutDistances <- distStitch(input$casType, input$cutSite)
						
					} else { # If no pre-made PAMs
						pams         <- pamStitch( "", input$customPamSeq)
						cutDistances <- distStitch("", input$cutSite)
					}
					)
				} else { # If no custom pams
					pams         <- input$casType
					cutDistances <- rep(-3, length(input$casType))
				}
				
				#if(input$talenOp == 1){
				#	talFlag <- TRUE
				#} else {
				talFlag <- FALSE
				#}
				#Calculate the MENTHU score
				#results <<- calculateMENTHUGeneSeqGenBank(pams, cutDistances, wiggle = TRUE, wiggleRoom = 39, talenList, gbFlag, gbhFlag, talFlag,
				#																					info, input$threshold, input$firstExon, input$exonTargetType, exonStuff, 
				#																					progress, input$scoreScheme)
				#results <<- calculateMENTHUGeneSeqGenBank(pams, cutDistances, wiggle = TRUE, wiggleRoom = 39, talenList, gbFlag, gbhFlag, 
				#																					info, input$threshold, input$firstExon, input$exonTargetType, exonStuff, progress)
				stuff    <- calculateMENTHUEnsembl(pams, cutDistances, wiggle = TRUE, wiggleRoom = 39, talenList, ensemblInfo, exonStuff, progress)
				
				
				# Statistics on the number of targets detected
				output$ensemblHits <- renderUI({
					HTML(paste(
						paste0("Number of target sites detected: ",                                             stuff[[2]]),
						paste0("Number of target sites with sufficient sequence context to calculate score: ",  stuff[[3]]),
						paste0("Number of target sites with 3bp microhomology arms within 5bp of each other: ", stuff[[4]]),
						paste0("Number of target sites with score above 1.5 threshold: ",                       stuff[[5]]),
						paste0("Number of target sites satisfying 3bp mh and threshold constraints: ",          stuff[[6]]),
						sep = "<br>"
					))
				})
				
				results <<- stuff[[1]]
				
				# Order the result table from largest menthuScore to smallest, and drop 0s
				results <<- results[which( results$MENTHU_Score > 0), ]
				results <<- results[order(-results$MENTHU_Score),     ]
				
				# Set the flag to display genbank results
				rValues$ensemblResultsFlag <- TRUE
				
				# Set the download button flag to true to render download button visible
				rValues$downloadFEns <- TRUE 
				print(paste0("Time to calculate: ", Sys.time() - timeX))
				
				} else {
					
				}
							
			} else {
				# Output error if no genbank file is found
				output$ensemblIdOutcome <- renderText(paste0("Error: Accession '", input$ensemblId, "' was not found in database."))
			}
			
		} else {
			# Output error if Ensembl is down
			output$ensemblUp <- renderText(paste0("Error: MENTHU was unable to connect to Ensembl."))
		}
		
	})
	
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
	#Copy/paste Submit Handler
	observeEvent(input$geneSeqSubmit, {
		
		#Run checks for okay PAM/TALEN input
		#talFlag <- 2
		#if(input$talenOp == 1){
		#	if(is.null(validTalen())){
		#		talFlag <- 1 #Input TALEN is good
		#	} else {
		#		talFlag <- 2 #Problems with input TALEN
		#	}
		#} else {
			talFlag <- 0 #TALENs not selected
		#}
		
		# Run checks for okay custom PAM input
		cusPamFlag <- 2
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) && is.null(validCustomCutSites())){
				cusPamFlag <- 1
				
			} else {
				cusPamFlag <- 2
				
			}
		} else {
			cusPamFlag <- 0
		}
		
		# Run checks for okay length match b/w PAM and cut site
		lenMatch <- 2
		if(input$customCutOpt == 1){
			if(is.null(validMatchCustomInputLength())){
				lenMatch <- 1
				
			} else {
				lenMatch <- 2
				
			}
		} else {
			lenMatch <- 0
		}
		
		# Prevent from running without okay inputs
		if(is.null(validGeneSeq()) &&           # Check Genbank ID is okay
			 #is.null(validThreshold()) &&         # Check threshold is okay
			 is.null(validPAM()) &&               # Check that one of the input options is selected
			 (cusPamFlag != 2) &&                 # Check if custom PAMs are used, and if they are okay
			 (talFlag != 2) &&                    # Check if TALENs are used, and if they are okay
			 (lenMatch != 2)){                     
			
			# Make a progress object
			progress <- shiny::Progress$new()
			
			# Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			# If there is not exon input, just set to 0, otherwise use exon input
			if(input$pasteExonType == 0){
				exonIn <- 0
				
			} else {
				if(is.null(input$exonInfo)){
					exonIn <- 0
					
				} else {
					exonIn <- exonHandler(input$exonInfo)
				}
			}
			
			# Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			if(talFlag == 0){
				talArmin  <- ""
				talArmax  <- ""
				talSpamin <- ""
				talSpamax <- ""
				
			} else {
				
				#if(input$talenOp == 1){
				#	if(input$spacer == 0){
				#		spamin <- 14
				#		spamax <- 14
				#	} else if(input$spacer == 1){
				#		spamin <- 16
				#		spamax <- 16
			#			} else {
				#			spamin <- 14
				#			spamax <- 16
				#		}
				#}
				#talArmin  <- input$armin
				#talArmax  <- input$armax
				#talSpamin <- spamin
				#talSpamax <- spamax
			}
			
			# Format custom PAMs and cut distances for submission to calculateMENTHUGeneSeq
			suppressWarnings(if(input$customCutOpt == 1){ # If customs pams in use
				if(is.null(input$casType)){
					tCasFlag <- FALSE
				} else {
					tCasFlag <- TRUE
				}
				
				if(tCasFlag){
					pams         <- pamStitch( input$casType, input$customPamSeq)
					cutDistances <- distStitch(input$casType, input$cutSite)
					
				} else { # If no pre-made PAMs
					pams         <- pamStitch( "", input$customPamSeq)
					cutDistances <- distStitch("", input$cutSite)
					
				}
				
			} else { # If no custom pams
				pams <- input$casType
				cutDistances <- rep(-3, length(input$casType))
			})
			
			#if(input$talenOp == 1){
			#	talFlag <- TRUE
			#} else {
				talFlag <- FALSE
			#}
			
			# Calculate the MENTHU score
			#results <<- calculateMENTHUGeneSeq(pams, cutDistances, wiggle = TRUE, wigRoom = 39, toupper(input$geneSeq), input$threshold, talFlag,
			#																	 exonIn, progress, talArmin, talArmax, talSpamin, talSpamax, input$scoreScheme)
			#results <<- calculateMENTHUGeneSeq(pams, cutDistances, wiggle = TRUE, wigRoom = 39, input$geneSeq, input$threshold, 
			#																	 exonIn, progress, talArmin, talArmax, talSpamin, talSpamax)
			
			stuff   <<- calculateMENTHUGeneSeq(pams, cutDistances, wiggle = TRUE, wigRoom = 39, input$geneSeq, exonIn, progress)
			
			results <<- stuff[[1]]
			
			if(is.numeric(results)){
				rValues$geneSeqError <- 1
			} else {
				
				# Output statistics regarding number of target sites found
				output$geneseqhits <- renderUI({
					HTML(paste(
						paste0("Number of target sites detected: ",                                             stuff[[2]]),
						paste0("Number of target sites with sufficient sequence context to calculate score: ",  stuff[[3]]),
						paste0("Number of target sites with 3bp microhomology arms within 5bp of each other: ", stuff[[4]]),
						paste0("Number of target sites with score above 1.5 threshold: ",                       stuff[[5]]),
						paste0("Number of target sites satisfying 3bp mh and threshold constraints: ",          stuff[[6]]),
						sep = "<br>"
					))
				})
				
				# Order the result table from largest menthuScore to smallest
				results <<- results[which(results$MENTHU_Score > 0), ]
				results <<- results[order(-results$MENTHU_Score),]
				
				rValues$geneSeqResultsFlag <- TRUE
				
				# Set the download button flag to true to render download button visible
				rValues$downloadF <- TRUE
			}
		}
	})
	
	########################################################
	#############Action Links/Reset Buttons#################
	########################################################
	
	####Update checkbox with select all####
	observeEvent(input$selectAll, {
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NGG", "NRG", "NNNRRT", "NNGRRT", "NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
	})
	
	####Update checkbox with select none (leaves SpCas9 NGG selected)
	observeEvent(input$selectNone, {
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = ""
		)
	})
	
	# Action link for inputting an example GenBank accession/session
	observeEvent(input$exampleGenBank, {
		reset()        # Clear all the inputs
		resetOutputs() # Clear the output regions
		
		# Reset/clear the input exon table
		rValues$rhFrame   <- dfEmpty
		rValues$resetVal  <- TRUE
		
		# Do not display the input exon table
		geneSeqExampleFlag <<- FALSE
		
		# Update the inputType to GenBank
		updateRadioButtons(session,
											 "inputType", 
											 selected = 1)
		
		# Update the PAM selection
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NRG", "NNNRRT"))
		
		# Pre-populate with the flh zebrafish gene
		updateTextAreaInput(session,
												"genbankId",
												value = "AY214391.1")
		
		# Update to use the first exon
		updateRadioButtons(session,
											 "firstExon",
											 selected = 1)
		
		# Search all exons
		updateRadioButtons(session,
											 "exonTargetType",
											 selected = 0)
		
		#Update the threshold if scoring scheme = version 1
		#if(input$scoreScheme == 1){
		#	updateNumericInput(session,
		#										 "threshold",
		#										 value = 40)
		#} else {
		#Update the threshold if scoring scheme = version 2
		#	updateNumericInput(session,
		#										 "threshold",
		#										 value = 1.5)
		#}
		
		# Use custom PAM input
		updateRadioButtons(session,
											 "customCutOpt",
											 selected = 1)
		
		# Make custom PAMs
		updateTextAreaInput(session,
												"customPamSeq",
												value = "NNNGCT")
		
		# Make custom PAM distance
		updateTextAreaInput(session,
												"cutSite",
												value = "-4")
	})
	
	# Example of gene sequence copy/paste input
	observeEvent(input$exampleGeneSeq, {
		reset()
		resetOutputs()
		
		rValues$rhFrame   <- dfExample
		rValues$resetVal  <- TRUE
			
		# Set flag to display 'exon' information
		geneSeqExampleFlag <<- TRUE
		
		# Change input type to copy/paste
		updateRadioButtons(session,
											 "inputType",
											 selected = 2)
		
		# Change to custom exon input
		updateRadioButtons(session,
											 "pasteExonType",
											 selected = 1)
		
		# Sample copy/paste gene input
		updateTextAreaInput(session, "geneSeq",
												value = paste0("CTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGAT",
																			 "ACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGC",
																			 "GCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGA",
																			 "CAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATACGCGTACCGCTAGC",
																			 "CAGGAAGAGTTTGTAGAAACGCAAAAAGGCCATCCGTCAGGATGGCCTTCTGCTTAGTTTG",
																			 "ATGCCTGGCAGTTTATGGCGGGCGTCCTGCCCGCCACCCTCCGGGCCGTTGCTTCACAACG",
																			 "TTCAAATCCGCTCCCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGA",
																			 "TAAAACGAAAGGCCCAGTCTTCCGACTGAGCCTTTCGTTTTATTTGATGCCTGGCAGTTCC",
																			 "CTACTCTCGCGTTAACGCTAGCATGGATGTTTTCCCAGTCACGACGT"))
		
		# Sample PAMs
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
		
		# Use custom PAM input
		updateRadioButtons(session,
											 "customCutOpt",
											 selected = 1)
		
		# Make custom PAMs
		updateTextAreaInput(session,
												"customPamSeq",
												value = "NNNGAG NNNGNG")
		
		# Make custom PAM distance
		updateTextAreaInput(session,
												"cutSite",
												value = "-4 -3")
		
		#if(input$scoreScheme == 1){
			#Change threshold to 30 to show output
		#	updateNumericInput(session, "threshold", value = 30)
		#} else {
			#Change threshold to 1.5 to show output
		#	updateNumericInput(session, "threshold", value = 1.5)
		#}
	})
	
	# When hitting the reset link button, call the reset function
	observeEvent(input$reset, {
		reset()
		resetOutputs()
	})
	
	resetOutputs <- function(){
		if(input$inputType == 1){
			resetGenBankOutputs()
		} else if(input$inputType == 2){
			resetGeneSeqOutputs()
		}
	}
	
	# Reset function
	reset <- function(){
		# Reset the geneSeqExample flag input
		geneSeqExampleFlag <<- FALSE
		
		# Update scoring system
		#updateRadioButtons(session,
		#									 "scoreScheme",
		#									 selected = 2)
		
		# Update checkboxes to only have SpCas9 NGG selected
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = "NGG"
		)
		
		# Clear custom pam input
		updateTextAreaInput(session,
												"customPamSeq",
												value = ""
		)
		
		# Clear custom DSB
		updateTextAreaInput(session,
												"cutSite",
												value = "")
		
		# Update custom PAM sequence to 'no'
		updateRadioButtons(session,
											 "customCutOpt",
											 selected = 0)

		if(input$inputType == 2){
			# Reset geneSeq input
			updateTextAreaInput(session, 
													"geneSeq", 
													label = "", 
													value = "")
			
			updateRadioButtons(session,
												 "pasteExonType",
												 selected = 0)

			# Reset genbank input
			#updateTextAreaInput(session,
			#										"genbankId",
			#										label = "",
			#										value = "")
			
		} else if(input$inputType == 1){
			# Reset genbank input
			updateTextAreaInput(session,
													"genbankId",
													label = "",
													value = "")
			#Reset geneSeq input
			#updateTextAreaInput(session, 
			#										"geneSeq", 
			#										label = "", 
			#										value = "")
			
			#updateRadioButtons(session,
			#									 "pasteExonType",
			#									 selected = 0)
<<<<<<< HEAD
=======
			
		} else if(input$inputType == 3){
			updateTextAreaInput(session,
													"ensemblId",
													label = "",
													value = "")
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
		}

		# Reset exon targets table
		if(input$inputType == 1 | input$inputType == 3){
			updateRadioButtons( session, "firstExon",         selected = 0)
			updateRadioButtons( session, "exonTargetType",    selected = 0)
			updateNumericInput( session, "exonBegPercentage", value    = 30)
			updateNumericInput( session, "exonEndPercentage", value    = 30)
			updateTextAreaInput(session, "exonTargetList",    value    = "")
		}
		
		
		# Reset sequence input type
		updateRadioButtons(session, 
											 "inputType", 
											 selected = 1
		)
		

		# Reset talen inputs
		#updateNumericInput(session, "spamin", value = 14)
		#updateNumericInput(session, "spamax", value = 16)
		#updateRadioButtons(session, "spacer", selected = 0)
		#updateNumericInput(session, "armin",  value = 15)
		#updateNumericInput(session, "armax",  value = 15)
		
		# Reset talen select
		#updateRadioButtons(session, "talenOp", selected = 0)
		
		# Reset threshold input
		#updateNumericInput(session, "threshold", value = 1.5)
		
<<<<<<< HEAD
=======
		# Reset intronic controls
		updateRadioButtons(session, "contextWiggleType", selected = 0)
		updateRadioButtons(session, "gRNAWiggleType",    selected = 0)
		
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
		# Reset rhandsontable
		rValues$resetVal <- TRUE
		rValues$rhFrame <- dfEmpty
		
		#
		rValues$validExonFlag <- TRUE
		
		# Set the download flag to not display the download button
<<<<<<< HEAD
		rValues$downloadF   <- FALSE
		rValues$downloadFGB <- FALSE
		
=======
		rValues$downloadF    <- FALSE
		rValues$downloadFGB  <- FALSE
		rValues$downloadFEns <- FALSE
	
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
		# Empty the result storage
		results <<- 0
	}
	
	resetGenBankOutputs <- function(){
		# Clear copy/paste outputs
		output$genbankIdOutcome <- renderText({
			""
		}) 
<<<<<<< HEAD
		
		rValues$genbankResultsFlag <- FALSE
		rValues$geneSeqResultsFlag <- FALSE	
=======
		clearOutputFlags()
	}
	
	# Clear the Ensembl unique outputs and also all output flags
	resetEnsemblOutputs <- function(){
		clearOutputFlags()
>>>>>>> parent of 7157e4f... Revert "Added Ensembl validation checks; updated download file name to include Ensembl/Genbank ID; updated resets to include Ensembl inputs"
	}
	
	resetGeneSeqOutputs <- function(){
		rValues$geneSeqResultsFlag <- FALSE	
		rValues$genbankResultsFlag <- FALSE
		rValues$geneSeqError <- 0
		output$genbankIdOutcome <- renderText({
			""
		})
	}
	
	# Functions for updating slider input based on user actions
	observeEvent(input$armin, {
		updateSliderInput(session = session, inputId = "armax",  min = input$armin)
	})
	
	observeEvent(input$armax, {
		updateSliderInput(session = session, inputId = "armin",  max = input$armax)
	})
	
	#observeEvent(input$spamin, {
	#	updateSliderInput(session = session, inputId = "spamax", min = input$spamin)
	#})
	

	#observeEvent(input$spamax, {
	#	updateSliderInput(session = session, inputId = "spamin", max = input$spamax)
	#})
})

