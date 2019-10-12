# Required packages
library(shiny)
library(shinyjs)
library(rhandsontable)
library(stringr)
library(stringi)
library(Biostrings)
library(rentrez)
library(rlist)
library(DT)
library(plyr)
library(Rcpp)
library(curl)
library(httr)
library(jsonlite)
library(xml2)

# Required supporting files
source("apeShiftFunctions.R")
source("genbankAccessoryFunctions.R")
source("menthu2.0AccessoryFunctions.R")
source("required2.0Functions_1.R")
source("targetAccessoryFunctions2.0.R")
source("ensemblAccessoryFunctions.R")

shinyServer(function(input, output, session){
	
	########################################################################################################################################################################
	###################Global Variables#####################################################################################################################################
	########################################################################################################################################################################
	
	# Storage for final results so that it is accessible to download handler and submit buttons
	results   <<- 0
	
	# Empty data frame for displaying when users are inputting exon information
	dfEmpty   <<- data.frame(Exon_Num         = rep(0, 5),
													 exonStart        = rep(0, 5), 
													 exonEnd          = rep(0, 5), 
													 stringsAsFactors = FALSE)
	
	pamEmpty   <<- data.frame(PAM_Sequence     = rep("NNN", 1),
													  DSB_Position     = rep(0, 1), 
													  Overhang_Length  = rep(0, 1), 
													  stringsAsFactors = FALSE)
	
	# Example data frame for when the paste gene seq example is selected
	dfExample <<- data.frame(Exon_Num         = c(1,     2,   3,   4,   5),
													 exonStart        = c(1,   101, 201, 301, 401),
													 exonEnd          = c(100, 200, 300, 400, 500),
													 stringsAsFactors = FALSE)
	
	# Example data frame for when examples are selected
	pamExample   <<- data.frame(PAM_Sequence     = c("NG", "TTN"),
													   	DSB_Position     = c(-3,   18), 
												  		Overhang_Length  = c(0,    5), 
												  		stringsAsFactors = FALSE)
	
	# Reactive values
	rValues <- reactiveValues(downloadF          = FALSE,    # Flag for displaying download button (for copy/paste gene seq)
														downloadFGB        = FALSE,    # Flag for displaying download button (for genbank)
														filtOpsGB          = FALSE,    # Flag for displaying filter options after calculation completion
														filtOpsGS          = FALSE,    # Flag for displaying filter options after calculation completion
														filtOpsE           = FALSE,    # Flag for displaying filter options after calculation completion
														geneSeqResultsFlag = FALSE,    # Flag for displaying gene seq results table
														genbankResultsFlag = FALSE,    # Flag for displaying genbank results table
														rhFrame            = dfEmpty,  # Slot to hold exon information data frame
														pamFrame           = pamEmpty, # Slot to hold custom nuclease information
														resetVal           = FALSE,    # For if the reset button has been clicked
														geneSeqError       = 0,        # Error messages for geneSeq submission
														downloadFEns       = FALSE,    # Flag for displaying download button (for Ensembl)
														ensemblResultsFlag = FALSE,    # Flag for displaying Ensembl results table
														validExonListFlag  = TRUE
	)
	
	# Load Ensembl ID table
	ensIds <<- readRDS("2018-09-21_ensIds.RDS")
	
	# Flag for displaying example table when clicking example geneSeq link
	geneSeqExampleFlag <<- FALSE
	
	########################################################################################################################################################################
	##################Validation Checks#####################################################################################################################################
	########################################################################################################################################################################
	
	####Make sure GenBank/RefSeq ID is properly formatted####
	validRefSeqGBGenbankId <- reactive({
		if(input$genbankId != ""){
			#Let RefSeq RNA accessions through
			if(       stringr::str_detect(input$genbankId, regex("^(NM|NR|XM|XR)_[0-9]{6}",       ignore_case = TRUE))){
				
				#Catch RefSeq protein accesssions
			} else if(stringr::str_detect(input$genbankId, regex("^(AP|NP|YP|XP|WP)_[0-9]{6}",    ignore_case = TRUE))){
				shiny::validate(
					need(1 == 2,
							 paste0("Error: This ID matches the RefSeq protein accession format. ", 
							 			 "Please submit an accession corresponding to a DNA sequence.")))
				
				# Catch ginormous genomic region files; disabled for running locally
			} else if(stringr::str_detect(input$genbankId, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}", ignore_case = TRUE))){
				#shiny::validate(
				#	need(1 == 2,
				#			 paste0("Error: This ID matches a RefSeq complete genomic molecule, incomplete genomic region, contig, ", 
				#			 			 "scaffold, or complete genome. We do not currently support any of these reference types due to ", 
				#			 			 "issues surrounding exon identification. You can use a GenBank or RefSeq nucleotide entry ", 
				#			 			 "corresponding to your gene of interest."))
				#)
				
				# Catch GenBank protein accessions
			} else if(stringr::str_detect(input$genbankId, regex("^[A-Z]{3}[0-9]{5}",             ignore_case = TRUE))){
				shiny::validate(
					need(1 == 2,
							 paste0("Error: This ID matches GenBank protein accession format. ", 
							 			 "Please submit an accession corresponding to a DNA sequence.")))
				
				# Catch anything else not conforming to input types
			} else {
				shiny::validate(
					need(((( stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{2}[0-9]{6}",                  ignore_case = TRUE)))  | 
								 	(stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{1}[0-9]{5}",                  ignore_case = TRUE)))) |
									(stringr::str_detect(input$genbankId, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}\\.[0-9]", ignore_case = TRUE)))) |
							 	(  stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{4}[0-9]{8,10}",               ignore_case = TRUE))),
							 paste0("Error: The entered ID does not match a known Genbank NUCLEOTIDE or RefSeq NUCLEOTIDE ID format.", 
							 			 "Please check your submitted accession to make sure you are using a NUCLEOTIDE entry.")))
			}
		} 
	})
	
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
	# If we can make contact with Ensembl
	ensemblIdExists <- reactive({
		if((input$inputType == 3) & input$ensemblId != ""){
			if(is.null(validEnsemblId())){
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
	
	####Validate copy/paste sequence input####
	validGeneSeq <- reactive({
		#If the input type is copy/paste gene seq and text has been entered
		if((input$geneSeq   != "") && 
			 (input$inputType == 2)){
			shiny::validate(
				#Check for fasta format
				need(!stringr::str_detect(input$geneSeq, "[>]"), 
						 "Error: Input DNA sequence appears to be in FASTA format. Please paste your sequence without the fasta header."),
				
				#Check for not DNA input
				need(!stringr::str_detect(input$geneSeq, "[^ACGTacgt0-9\\s\\n]"), 
						 paste0("Error: Input DNA sequence contains non-standard nucleotides. ", 
						 			 "Allowed nucleotides are A, C, G, and T.")),
				
				# Prevent users from submitting too short sequences
				need(nchar(input$geneSeq) >= 80,
						 paste0("The DNA sequence has <80 nucleotides. MENTHU requires at least 40 nucleotides upstream", 
						 			 "and 40 nucleotides downstream of the DSB site in order to properly calculate the MENTHU score."))#,
				
				#Prevent users from blowing up the server
				#need(nchar(input$geneSeq) < 5000, 
				#		 paste0("The DNA sequence has >5000 nucleotides. For sequences of this size, ",
				#		 			  "please use the local version of MENTHU, which can be accessed via the 'Tools and Downloads' tab."))
			)
			
		} else if(input$inputType == 2) {
			# Prevent running on empty submission
			shiny::validate(
				need(input$geneSeq != "", "")
			)
		}
	})
	
	# Valid PAM - make sure user selects at least one target type
	validPAM <- reactive({
		shiny::validate(
			need((length(input$casType) > 0) | (input$talenOp == 1) | (input$customCutOpt == 1),
					 "Error: No nuclease selected. Please select a Cas type and/or a custom PAM scheme and/or the TALEN input option.")
		)
	})
	
	# Valid TALEN inputs
	validTalen <- reactive({
		shiny::validate(
			# Make sure that the arm length min is less than max
			need(input$armin   <= input$armax,
					 "Error: Maximum TALEN arm length must be greater than or equal to minimum TALEN arm length."),
			# Limit arm length to 15-18 nt
			need((input$armin  >= 15),
					 "Error: Minimum TALEN arm length is 15 nucleotides."),
			need((input$armax  <= 18),
					 "Error: Maximum TALEN arm length is 18 nucleotides.")
		)
	})
	
	# Valid custom PAM inputs
	validCustomPam <- reactive({
		# If the custom PAM is selected, make sure it only has IUPAC nucleotides or separating characters
		if(input$customPamSeq != ""){
			if(stringr::str_detect(input$customPamSeq, "[^ACGTRYSWKMBDHVNacgtryswkmbdhvn\\s,]")){
				shiny::validate(
					need(1 == 2,
							 paste0("Error: Non-allowed characters detected. ", 
							 			  "Only standard nucleotide (A, C, G, T), IUPAC extended ", 
							 			  "nucleotide symbols (R, Y, S, W, K, M, B, D, H, V, N), ", 
							 			  "and separating characters (spaces and commas) allowed."))
				)
			} else {
				# Check for all 'N' PAMs, single nucleotide PAMs, or PAMs consisting solely of Ns and a single nucleotide; disabled for running locally
				#disallowed      <- "^[N]+[ACGTRYSWKMBDHV]{1}$|^[N]+[ACGTRYSWKMBDHV]{1}[N]+$|^[ACGTRYSWKMBDHVN]{1}$|^[N]+$"
				#disallowedCheck <- grepl(disallowed, pamStitch("", input$customPamSeq), ignore.case = TRUE, perl = TRUE)
				#disallowedExist <- is.element(TRUE, unlist(disallowedCheck))
				
				#if(disallowedExist){
				#	shiny::validate(
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
			shiny::validate(
				need(!stringr::str_detect(input$cutSite, "[^0-9\\s,-]"),
						 "Error: Only numbers, dashes, commas, and spaces allowed.")
			)
		}
	})
	
	# Valid custom overhangs
	validOverhangs <- reactive({
		if(input$overhang != ""){
			shiny::validate(
				need(!stringr::str_detect(input$overhang, "[^0-9\\s,-]"),
						 "Error: Only numbers, dashes, commas, and spaces allowed.")
			)
		}
	})
	
	# Check that the number of custom PAMs matches the number of custom cutSites
	validMatchCustomInputLength <- reactive({
		if((input$customPamSeq != "") && 
			 (input$customCutOpt != "")){
			shiny::validate(
				need(length(as.character(distStitch("", input$cutSite))) == length(pamStitch("", input$customPamSeq)) &&
						 length(as.character(distStitch("", input$cutSite))) == length(as.character(ohStitch("", input$overhang))),
						 paste0("Error: The number of custom PAM sequences, custom DSB sites, and overhang lengths does not match. ", 
						 			  "Please make sure that each PAM has one distance to its cut site and overhang length specified."))
			)
		}
	})
	
	# Check to make sure that exon list inputs are interpretable
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
	output$validgenbankid              <- renderText({validRefSeqGBGenbankId()     })
	
	output$genbankidexists             <- renderText({                             })
	
	output$validensemblid              <- renderText({validEnsemblId()             })
	
	output$ensemblidexists             <- renderText({ensemblIdExists()            })
	
	output$validgeneseq                <- renderText({validGeneSeq()               })
	
	output$validpam                    <- renderText({validPAM()                   })
	
	output$validtalen                  <- renderText({validTalen()                 })
	
	output$validcustompam              <- renderText({validCustomPam()             })
	
	output$validcustomcutsites         <- renderText({validCustomCutSites()        })
	
	output$validoverhangs              <- renderText({validOverhangs()             })
	
	output$validmatchcustominputlength <- renderText({validMatchCustomInputLength()})
	
	output$geneseqerrors               <- renderText({geneSeqErrors()              })
	
	output$validexonlist               <- renderText({validExonList()              })
	
	#output$validexoninfo               <- renderText({validExonInfo()              })
	
	########################################################################################################################################################################
	##################Bookmark Functions####################################################################################################################################
	########################################################################################################################################################################
	
	observeEvent(input$bookmarkGS, {session$doBookmark()})
	
	observeEvent(input$bookmarkGB, {session$doBookmark()})
	
	observeEvent(input$bookmarkE,  {session$doBookmark()})
	
	########################################################################################################################################################################
	#####################UI Rendering#######################################################################################################################################
	########################################################################################################################################################################
	
	# If there are no PAMs in the target sequence
	geneSeqErrors <- reactive({
		rValues$geneSeqError
		if(rValues$geneSeqError == 0){
			""
		} else if(rValues$geneSeqError == 1){
			"Error: No targets corresponding to the PAM(s) detected."
		}
	})
	
	# TODO - do we need this function for TALENs?
	
	
	
	###########################################################
	######### Rendering functions for exon input table ########
  #########           ROLL FOR INITIATIVE         ###########
	###########################################################
	#                     ,     \    /      ,                 #
  #	                   / \    )\__/(     / \                #
	#						        /   \  (_\  /_)   /   \               #
	#						   ____/_____\__\@  @/___/_____\____          #
	#					  	|             |\../|              |         #
	#						 	|              \VV/               |         #
	#						 	|          Here thar be           |         #
	#             |            DARGONS!!            |         #
	#						 	|_________________________________|         #
	#						 	|    /\ /      \\       \ /\    |           #
	#						 	|  /   V        ))       V   \  |           #
	#             |/     `       //        '     \|           #
  #             `              V                '           #
	###########################################################
	# Dragon ASCII art from                                   #
	# https://www.asciiart.eu/mythology/dragons               #
	###########################################################
	# But for realsies do not touch any of this; it literally #
	# took months to get lines 397-435 working properly.      #
	###########################################################
	
	#### Render function for exon input table
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
	
	#### Render exon input table
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
	
	#### Render table for nucleases ####
	renderTableNuc <- reactive({
		pamOutDF <- NULL
		
		if(!is.null(input$pamTable)){
			pamOutDF <- hot_to_r(input$pamTable)
			
		} else if(!is.null(isolate(rValues$pamFrame))){
			pamOutDF <- isolate(rValues$pamFrame)
		}
		
		if(!is.null(pamOutDF)){
			
			rValues$pamFrame <- pamOutDF
		}
		
		return(pamOutDF)
	}) %>% debounce(1000)
	
	#### Render more nuclease table ####
	output$pamTable <- renderRHandsontable({
		# Action buttons that this function is dependent on
		input$reset
		input$exampleGeneSeq
		input$exampleGenBank
		input$exampleEnsembl
		
		if(isolate(rValues$resetVal)){
			pamOutDF <- rValues$pamFrame
		} else {
			pamOutDF <- renderTableNuc()
			rValues$resetVal <- FALSE
		}
		
		if(!is.null(pamOutDF)){
			rhandsontable(pamOutDF) %>%
				hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>% 
				hot_col("PAM_Sequence",    format = "N") %>%
				hot_col("DSB_Position",    format = 0) %>%
				hot_col("Overhang_Length", format = 0)
		}
	})
	
	##########################################################
	#################### END DARGONS #########################
	##########################################################
	
	##########################################################
	################# Download Stuff #########################
	##########################################################
	
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
	
	# Download button for Ensembl results
	output$downOutEns <- renderUI({
		if(rValues$downloadFEns){
			downloadButton("downRes", "Download Results")
		} else {
			""
		}
	})
	
	# Download handler code
	output$downRes <- downloadHandler(
		filename = function(){
			#Name file in the form "YYYY-MM-DD_HH-MM-SS_InputID_targets.csv
			if(input$inputType == 1){
				# For GenBank/Refseq
				# Tag the file if the outputs were filtered
				if(input$t7OptGB || input$thresholdGB || !input$inFrameGB || !input$outFrameGB){
					filt <- "filtered"
				} else {
					filt <- ""
				}
				# Generate the file name
				paste0(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_", input$genbankId, "_MENTHU_targets_", filt, ".csv")
				
			} else if(input$inputType == 3){
				# For Ensembl
				# Tag the file if the outputs were filtered
				if(input$t7OptE || input$thresholdE || !input$inFrameE || !input$outFrameE){
					filt <- "filtered"
				} else {
					filt <- ""
				}
				# Generate the file name
				paste0(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_", input$ensemblId, "_MENTHU_targets_", filt, ".csv")
				
			} else {
				# For copy/paste
				# Tag the file if the outputs were filtered
				if(input$t7OptGS || input$thresholdGS || !input$inFrameGS || !input$outFrameGS){
					filt <- "filtered"
				} else {
					filt <- ""
				}
				# Generate the file name
				paste0(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_custom_seq_MENTHU_targets_", filt, ".csv")
				
			}
		},
		
		content = function(file){
			# Apply filters to the download
			if(input$inputType == 1){
				resOut <- filterResults(results, input$t7OptGB, input$thresholdGB, input$inFrameGB, input$outFrameGB)
				
			} else if(input$inputType == 2){
				resOut <- filterResults(results, input$t7OptGS, input$thresholdGS, input$inFrameGS, input$outFrameGS)
				
			} else if(input$inputType == 3){
				resOut <- filterResults(results, input$t7OptE,  input$thresholdE, input$inFrameE, input$outFrameE)
				
			}
			
			# Check to make sure there's actually output to download
			
			# TODO allow users to remove PAMs from download
			if(resOut[[1]]){
				# Remove HTML comments
				resOut[[2]]$Target_Sequence <- gsub("<strong>|</strong>", "", resOut[[2]]$Target_Sequence, ignore.case = TRUE, perl = TRUE)
				
				# Output the file
				write.table(resOut[[2]], file, row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = ",")
			} else {
				
			}
		}
	)
	
	######################################################################################################################################################################
	################# UI FILTER OPTIONS ##################################################################################################################################
	######################################################################################################################################################################
	
	# Display filter conditional panel when output is finished calculating
	output$filtOpsGS <- reactive({return(rValues$filtOpsGS)})
	output$filtOpsGB <- reactive({return(rValues$filtOpsGB)})
  output$filtOpsE  <- reactive({return(rValues$filtOpsE )})
  
  # Make sure the filter options stay awake when hidden
  outputOptions(output, "filtOpsGS", suspendWhenHidden = FALSE)
  outputOptions(output, "filtOpsGB", suspendWhenHidden = FALSE)
  outputOptions(output, "filtOpsE",  suspendWhenHidden = FALSE)
	
	#####################################################################################################################################################################
	############################# RESULT OUTPUT #########################################################################################################################
	#####################################################################################################################################################################
	
  # Output function for copy/paste results
  output$geneSeqResults <- renderUI({
  	
  	rValues$geneSeqResultsFlag
  	input$inputType
  	
  	if(rValues$geneSeqResultsFlag && input$inputType == "2"){
  		out <- filterResults(results, input$t7OptGS, input$thresholdGS, input$inFrameGS, input$outFrameGS)
  		
  		if(out[[1]]){
  			output$placeholder <- DT::renderDataTable(out[[2]], 
  																								options  = list(scrollX = TRUE), 
  																								rownames = FALSE,
  																								escape   = FALSE)
  			DT::dataTableOutput("placeholder")	
  		} else {
  			"No sites satisfy the selected filters."
  		}
  		
  	} else {
  		""
  	}
  })
  
  # Output function for Genbank results
  output$genbankResults <- renderUI({
		rValues$genbankResultsFlag
  	input$inputType
  	
  	if(rValues$genbankResultsFlag && input$inputType == "1"){
  		out <- filterResults(results, input$t7OptGB, input$thresholdGB, input$inFrameGB, input$outFrameGB)
  		
  		if(out[[1]]){
  			output$placeholder <- DT::renderDataTable(out[[2]], 
  																								options  = list(scrollX = TRUE), 
  																								rownames = FALSE,
  																								escape   = FALSE)
  			DT::dataTableOutput("placeholder")	
  		} else {
  			"No sites satisfy the selected filters."
  		}
  		
  	} else {
  		""
  	}
  })
  
  # Output function for Ensembl results
  output$ensemblResults <- renderUI({

  	rValues$ensemblResultsFlag
  	input$inputType
  	
  	if(rValues$ensemblResultsFlag && input$inputType == "3"){
  		out <- filterResults(results, input$t7OptE, input$thresholdE, input$inFrameE, input$outFrameE)
  		
  		if(out[[1]]){
  			output$placeholder <- DT::renderDataTable(out[[2]], 
  																								options  = list(scrollX = TRUE), 
  																								rownames = FALSE,
  																								escape   = FALSE)
  			DT::dataTableOutput("placeholder")	
  		} else {
  			"No sites satisfy the selected filters."
  		}
  		
  	} else {
  		""
  	}
  })
	
	########################################################################################################################################################################
	#################Submission Handling####################################################################################################################################
	########################################################################################################################################################################
	
  ########################################################################################################################################################################
  ############### GenBank/RefSeq Code#####################################################################################################################################
  ########################################################################################################################################################################
	observeEvent(input$genbankSubmit,{
		resetOutputs()
		results <<- 0
		
		# Run checks for okay PAM/TALEN input
		if(input$talenOp == 1){
			if(is.null(validTalen())){
				talFlag <- 1 # TALEN input is good
			} else {
				talFlag <- 2 # Problems with input TALEN
			}
		} else {
			  talFlag <- 0 # TALENs not selected
		}
		
		# Check for valid custom PAMs
		# Set custom PAM to throw a problem unless corrected
		cusPamFlag <- 2
		
		# Determine custom PAM cut site validation flag
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) &&     
				 is.null(validCustomCutSites())){ 
				cusPamFlag <- 1 # Custom PAM cut sites are valid
				
			} else {
				cusPamFlag <- 2 # Custom PAM cut sites are NOT valid
				
			}
			
		} else {
			  cusPamFlag <- 0 # Custom PAM cut sites are not used
		}
		
		# Check for valid overhang flags
		# Set custom overhang to throw a flag unless corrected
		cusOhFlag <- 2
		
		# Determine the overhang validation falg
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) && 
				 is.null(validOverhangs())){
				cusOhFlag <- 1 # Custom overhangs are valid
				
			} else {
				cusOhFlag <- 2 # Custom overhangs are NOT valid
			}
		} else {
			  cusOhFlag <- 0 # Custom overhangs are not used
		}
		
		# Check to make sure that there are equal numbers of PAMs, cut sites, and overhangs
		lenMatch <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validMatchCustomInputLength())){
				lenMatch <- 1 # Correct number
				
			} else {
				lenMatch <- 2 # Incorrect number
				
			}
		} else {
			  lenMatch <- 0 # Custom PAMs aren't used
		}
		
		# Prevent the whole shebang from running without okay inputs
		if( is.null(validRefSeqGBGenbankId()) && # Check Genbank ID is okay
			  is.null(validPAM())               && # Check that one of the input options is selected
			 (cusPamFlag != 2)                  && # Check if custom PAMs are used, and if they are okay
			 (cusOhFlag  != 2)                  && # Check that the custom overhang flag is okay
			 (talFlag    != 2)                  && # Check if TALENs are used, and if they are okay
			 (lenMatch   != 2)){                   # Check the PAM/cutsite/overhang numbers 
			
			# Create a new progress object
			progress <- shiny::Progress$new()
			
			# Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			# Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			talenList <- ""
			
			# If talen options are used, get the spacer min/max values
			if(input$talenOp == 1){
				if(input$spacer == 0){
					spamin <- 14
					spamax <- 14
					
				} else if(input$spacer == 1){
					spamin <- 16
					spamax <- 16
					
				} else {
					spamin <- 14
					spamax <- 16
					
				}
				
				talenList <- list(input$armin, input$armax, spamin, spamax)
			}
			
			# Update progress
			progress$set(detail = "Retrieving GenBank entry...", value = 0.1)
			
			# Try to pull genbank entry associated with accession
			# Get the GenBank sequence with exon/intron information
			gba <- input$genbankId
			
			# TODO is all this still necessary?
			endOfTry <<- FALSE
			gbFlag   <<- FALSE
			gbhFlag  <<- FALSE
			
			#Try to retrieve the Genbank file 
			if(endOfTry == FALSE){
				
				tryCatch({
					#info <- suppressWarnings(wonkyGenBankHandler(gba))
					info <- suppressWarnings(getGenbankFile(gba))
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
					
				} else if(input$exonTargetType  == 2){
					exonStuff <- input$exonEndPercentage / 100
					
				} else if(input$exonTargetType  == 3){
					exonStuff <- input$exonTargetList	
				}
				
				if(!is.null(input$casType)){
					preGen <- getPreGenPamList(input$casType)
				}
				
				# Handle cut distances and PAMs for input to calculateMENTHUGeneSeqGenBank
				if(input$customCutOpt == 1){ # If customs pams in use
					suppressWarnings(if(!is.null(input$casType)){ # If pre-made PAMs in use
						pams         <- pamStitch( preGen$pamList, input$customPamSeq)
						cutDistances <- distStitch(preGen$cutDist, input$cutSite)
						overhangs    <- ohStitch(  preGen$ohList,  input$overhang)
						
					} else { # If no pre-made PAMs
						pams         <- pamStitch( "", input$customPamSeq)
						cutDistances <- distStitch("", input$cutSite)
						overhangs    <- ohStitch(  "", input$overhang)
					})
					
				} else if(length(input$casType) > 0) { # If no custom pams AND preGen used
					pams         <- preGen$pamList
					cutDistances <- preGen$cutDist
					overhangs    <- preGen$ohList
				
				# If only TALENs are used
				} else {
					pams         <- NULL
					cutDistances <- NULL
					overhangs    <- NULL
					
				}
				
				if(input$talenOp == 1){
					talFlag <- TRUE
					
				} else {
				  talFlag <- FALSE
				}
				
				#Calculate the MENTHU score
				stuff    <- calculateMENTHUGeneSeqGenBank(pams, cutDistances, overhangs, wiggle = TRUE, wiggleRoom = 39, talenList, gbFlag, gbhFlag, 
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
				
				# Set flags to display output-dependent UI elements
				rValues$filtOpsGB          <- TRUE # Show filtering options when there is output to filter
				rValues$genbankResultsFlag <- TRUE # Display the results table
				rValues$downloadFGB        <- TRUE # Make the GenBank download button visible
				
				# Order the result table from largest menthuScore to smallest, and drop 0s
				results <<- results[which(results$MENTHU_Score > 0), ]
				results <<- results[order(-results$MENTHU_Score)   , ]
				

				
			} else {
				# Output error if no genbank file is found
				output$genbankIdOutcome <- renderText(paste0("Error: Accession '", input$genbankId, "' was not found in database."))
			}
		}
	})
	
  ########################################################################################################################################################################
  ############### Ensembl Code############################################################################################################################################
  ########################################################################################################################################################################
  
	observeEvent(input$ensemblSubmit,{
		resetOutputs()
		results <<- 0
		#timeX <- Sys.time()
		
		# Set TALEN flag to throw error without correction
		talFlag <- 2
		
		#Run checks for okay PAM/TALEN input
		if(input$talenOp == 1){
			if(is.null(validTalen())){
				talFlag <- 1 # Input TALEN is valid
				
			} else {
				talFlag <- 2 # Input TALEN is NOT valid
				
			}
		} else {
		    talFlag <- 0 # TALENs not selected
		    
		}
		
		# Check that all the custom PAM stuff is okay
		
		cusPamFlag <- 2
		cusOhFlag  <- 2
		lenMatch   <- 2
		
		if(input$customCutOpt == 1){
			# Check that custom DSBs are OK
			# Flag as not valid unless fixed
			if(is.null(validCustomPam()) && 
				 is.null(validCustomCutSites())){
				cusPamFlag <- 1 # Custom DSBs are valid
				
			} else {
				cusPamFlag <- 2 # Custom DSBs are NOT valid
				
			}
			
			# Check that custom overhangs are OK
			# Flag as not valid unless fixed
			if(is.null(validCustomPam()) && 
				 is.null(validOverhangs())){
				cusOhFlag <- 1 # Custom overhangs are valid
				
			} else {
				cusOhFlag <- 2 # Custom overhangs are NOT valid
				
			}
			
			# Check that there are correct number of custom PAM sequences, DSBs, and overhangs
			# Flag as not valid unless fixed
			if(is.null(validMatchCustomInputLength())){
				lenMatch <- 1 # Custom PAM inputs are correctly matched
				
			} else {
				lenMatch <- 2 # Custom PAM inputs are NOT correctly matched
				
			}
			
		} else {
			cusPamFlag <- 0 # Custom DSBs are not used
			cusOhFlag  <- 0 # Custom overhangs are not used
			lenMatch   <- 0 # Custom PAMs are not used
		}
		
		# Prevent the whole shebang from running without okay inputs
		if(is.null(validEnsemblId())  && # Check Ensembl ID is okay
			 is.null(ensemblIdExists()) && # Make sure the ID actually exists
			 is.null(validPAM())        && # Check that one of the input options is selected
			 (cusPamFlag != 2)          && # Check if custom PAMs are used, and if they are okay
			 (cusOhFlag  != 2)          && # Check that the overhangs are ok
			 (talFlag    != 2)          && # Check if TALENs are used, and if they are okay
			 (lenMatch   != 2)){                    
			
			# Create a new progress object
			progress <- shiny::Progress$new()
			
			# Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			# Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			talenList <- ""
			
			# Check if TALENs are used, and set spacer min/max
			if(input$talenOp == 1){
				if(input$spacer == 0){
					spamin <- 14
					spamax <- 14
					
				} else if(input$spacer == 1){
					spamin <- 16
					spamax <- 16
					
				} else {
					spamin <- 14
					spamax <- 16
					
				}
				
				talenList <- list(input$armin, input$armax, spamin, spamax)
			}
			
			# Update progress when Ensembl retrieval starts
			progress$set(detail = "Retrieving Ensembl entry...", value = 0.1)
			
			# Make sure Ensembl is up and responsive
			if(isEnsemblUp()){
				ensemblInfo <- handleEnsemblInput(input$ensemblId, wiggle = TRUE, wiggleRoom = 39)
				
				progress$set(detail = "Processing Ensembl entry...", value = 0.1)
				
				# If the entry is NOT an exon
				if(getEnsemblIdType(input$ensemblId, check = TRUE) != "exon"){
					
					# Figure out how many exons there are in the transcript/protein
					ensemblInfo$rank <- as.numeric(ensemblInfo$rank)
					
					# Determine how many exons are present
					maxR <- max(ensemblInfo$rank)
					
					# Set the 'first' exon to be examined to exon 1 or exon 2 depending on factors
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
					
					if(!is.null(input$casType)){
						preGen <- getPreGenPamList(input$casType)
					}	
					
					# TODO Deal with Cas12a
					# Handle cut distances and PAMs for input to calculateMENTHUEnsembl
					if(input$customCutOpt == 1){ # If custom PAMs in use
						suppressWarnings(if(!is.null(input$casType)){ # If pre-made PAMs are also in use
							pams         <- pamStitch( preGen$pamList, input$customPamSeq)
							cutDistances <- distStitch(preGen$cutDist, input$cutSite)
							overhangs    <- ohStitch(  preGen$ohList,  input$overhang)
							
							
						} else { # If no pre-made PAMs
							pams         <- pamStitch( "", input$customPamSeq)
							cutDistances <- distStitch("", input$cutSite)
							overhangs    <- ohStitch(  "", input$overhang)
						})
						
					} else if(length(input$casType > 0)){ # If no custom pams
						pams         <- preGen$pamList
						cutDistances <- preGen$cutDist
						overhangs    <- preGen$ohList
					} else {
						pams         <- NULL
						cutDistances <- NULL
						overhangs    <- NULL
					}
					
					
					# Flag for whether or not TALENs are used
					if(input$talenOp == 1){
						talFlag <- TRUE
					} else {
					  talFlag <- FALSE
					}
					
					#Calculate the MENTHU score
					stuff    <- calculateMENTHUEnsembl(pams, cutDistances, overhangs, wiggle = TRUE, wiggleRoom = 39, talenList, ensemblInfo, exonStuff, progress)
					
					# Statistics on the number of targets detected
					output$ensemblHits <- renderUI({
						HTML(paste(
							paste0("Number of target sites detected: ",                                             stuff[[2]]),
							paste0("Number of target sites with sufficient sequence context to calculate score: ",  stuff[[3]]),
							paste0("Number of target sites with 3bp microhomology arms within 5bp of each other: ", stuff[[4]]),
							paste0("Number of target sites with score above 1.5 threshold: ",                       stuff[[5]]),
							paste0("Number of target sites satisfying 3bp mh and threshold constraints: ",          stuff[[6]]),
							sep = "<br />"
						))
					})
					
					results <<- stuff[[1]]
					
					# Set flags to display output-dependent UI elements
					rValues$filtOpsE           <- TRUE # Display filter options
					rValues$ensemblResultsFlag <- TRUE # Set the flag to display ensembl results
					rValues$downloadFEns       <- TRUE # Make the download button visible
					
					# Order the result table from largest menthuScore to smallest, and drop 0s
					results <<- results[which( results$MENTHU_Score > 0), ]
					results <<- results[order(-results$MENTHU_Score)    , ]
					
					#print(paste0("Time to calculate: ", Sys.time() - timeX))
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
	
  ########################################################################################################################################################################
  ############### Copy/Paste Handler Code ################################################################################################################################
  ########################################################################################################################################################################
  
	#Copy/paste Submit Handler
	observeEvent(input$geneSeqSubmit, {
		resetOutputs()
		results <<- 0
		# Run checks for okay PAM/TALEN input
		# Flag as not valid unless fixed
		talFlag <- 2
		
		if(input$talenOp == 1){
			if(is.null(validTalen())){
				talFlag <- 1 # TALEN inputs are valid
			} else {
				talFlag <- 2 # TALEN inputs are NOT valid
			}
		} else {
		    talFlag <- 0 # TALENs are not used
		}
		
		# Run checks for valid custom PAM DSBs
		# Flag as not valid unless fixed
		cusPamFlag <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) && 
				 is.null(validCustomCutSites())){
				cusPamFlag <- 1 # Custom PAM DSBs are valid
				
			} else {
				cusPamFlag <- 2 # Custom PAM DSBs are not valid
			}
			
		} else {
			  cusPamFlag <- 0 # Custom PAMs are not used
		}
		
		# Run checks for valid custom PAM overhangs
		# Flag as not valid unless fixed
		cusOhFlag <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validCustomPam()) && 
				 is.null(validOverhangs())){
				cusOhFlag <- 1 # Custom overhangs are valid
				
			} else {
				cusOhFlag <- 2 # Custom overhangs are not valid
			}
			
		} else {
		   	cusOhFlag <- 0 # Custom overhangs are not used
			
		}
		
		# Run checks for okay length match b/w PAM and cut site
		# Flag as not valid unless fixed
		lenMatch <- 2
		
		if(input$customCutOpt == 1){
			if(is.null(validMatchCustomInputLength())){
				lenMatch <- 1 # Correct number of PAM seqs, DSBs, and overhangs
				
			} else {
				lenMatch <- 2 # Invalid matchup
			}
			
		} else {
		  	lenMatch <- 0 # Customs PAMs are not used
		}
		
		# Prevent from running without okay inputs
		if(is.null(validGeneSeq()) &&           # Check Genbank ID is okay
			 is.null(validPAM())     &&           # Check that one of the input options is selected
			 (cusPamFlag != 2)       &&           # Check if custom PAMs are used, and if they are okay
			 (cusOhFlag  != 2)       &&           # Check if custom overhangs are valid
			 (talFlag    != 2)       &&           # Check if TALENs are used, and if they are okay
			 (lenMatch   != 2)){                     
			
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
				# Set the TALEN arm lengths and spacer lengths
				if(input$talenOp == 1){
					if(input$spacer == 0){
						spamin <- 14
						spamax <- 14
						
					} else if(input$spacer == 1){
						spamin <- 16
						spamax <- 16
						
					} else {
						spamin <- 14
						spamax <- 16
					}
				}
				
				talArmin  <- input$armin
				talArmax  <- input$armax
				talSpamin <- spamin
				talSpamax <- spamax
			}
			
			if(!is.null(input$casType)){
				preGen <- getPreGenPamList(input$casType)
			}
			
			# TODO Fix this for Cas12a
			# Format custom PAMs and cut distances for submission to calculateMENTHUGeneSeq
			suppressWarnings(if(input$customCutOpt == 1){
				
				# If premade PAMs are not used
				if(is.null(input$casType)){
					pams         <- pamStitch( "", input$customPamSeq)
					cutDistances <- distStitch("", input$cutSite)
					overhangs    <- ohStitch(  "", input$overhang)

					# Premade PAMs are used
				} else { 
					pams         <- pamStitch( preGen$pamList, input$customPamSeq)
					cutDistances <- distStitch(preGen$cutDist, input$cutSite)
					overhangs    <- ohStitch(  preGen$ohList,  input$overhang)
					
				}
				
				# No custom PAMs
			} else if(!is.null(input$casType)){ 
				pams         <- preGen$pamList
				cutDistances <- preGen$cutDist
				overhangs    <- preGen$ohList

			} else {
				pams         <- NULL
				cutDistances <- NULL
				overhangs    <- NULL
			}
			)
			
			if(input$talenOp == 1){
				talFlag <- TRUE
			} else {
			  talFlag <- FALSE
			}
			
			# Calculate the MENTHU score
			stuff   <<- calculateMENTHUGeneSeq(pams, cutDistances, overhangs, wiggle = TRUE, wiggleRoom = 39, 
																				 stripWhiteSpace(input$geneSeq), exonIn, progress, talArmin, talArmax, talSpamin, talSpamax)

			results <<- stuff[[1]]
			
			if(is.numeric(results)){
				rValues$geneSeqError <- 1
				
			} else {
				
				# Output statistics regarding number of target sites found
				# output$geneseqhits <- renderUI({
				# 	HTML(paste(
				# 		paste0("Number of target sites detected: ",                                             stuff[[2]]),
				# 		paste0("Number of target sites with sufficient sequence context to calculate score: ",  stuff[[3]]),
				# 		paste0("Number of target sites with 3bp microhomology arms within 5bp of each other: ", stuff[[4]]),
				# 		paste0("Number of target sites with score above 1.5 threshold: ",                       stuff[[5]]),
				# 		paste0("Number of target sites satisfying 3bp mh and threshold constraints: ",          stuff[[6]]),
				# 		sep = "<br>"
				# 	))
				# })
				
				# Show output-dependent UI elements
				rValues$filtOpsGS          <- TRUE # Show output filtering options
				rValues$geneSeqResultsFlag <- TRUE # Show the output table
				rValues$downloadF          <- TRUE # Make the output table visible
				
				# Order the result table from largest menthuScore to smallest
				results <<- results[which(results$MENTHU_Score > 0), ]
				results <<- results[order(-results$MENTHU_Score),]
				
			}
		}
	})
	
	########################################################################################################################################################################
	#############Action Links/Reset Buttons#################################################################################################################################
	########################################################################################################################################################################
	
	# Select All pre-gen PAMs
	observeEvent(input$selectAll, {
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NGG", "NRG", "NNNRRT", "NNGRRT", "NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT", "TTTN", "TTTV", "TTN", "YTN"))
	})
	
	# Select None pre-gen PAMs
	observeEvent(input$selectNone, {
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = ""
		)
	})
	
	####################################
	######## EXAMPLE: GenBank ##########
	####################################
	observeEvent(input$exampleGenBank, {
		reset()        # Clear all the inputs
		resetOutputs() # Clear the output regions
		
		# Reset/clear the input exon table
		rValues$rhFrame   <- dfEmpty
		rValues$pamFrame  <- pamExample
		rValues$resetVal  <- TRUE
		
		# Do not display the input exon table
		geneSeqExampleFlag <<- FALSE
		
		#### Update functions ####
		updateCheckboxGroupInput(session, "casType",        selected = c("NRG", "TTTN"))   # Update the PAM selection
		updateRadioButtons(      session, "inputType",      selected = 1)                  # Update the inputType to GenBank
		updateRadioButtons(      session, "firstExon",      selected = 1)                  # Update to use the first exon
		updateRadioButtons(      session, "exonTargetType", selected = 0)                  # Search all exons
		updateRadioButtons(      session, "customCutOpt",   selected = 1)                  # Use custom PAM input
		updateTextAreaInput(     session, "genbankId",      value    = "AY214391.1")       # Pre-populate with the flh zebrafish gene
		updateTextAreaInput(     session, "customPamSeq",   value    = "NG YYN")           # Make custom PAMs
		updateTextAreaInput(     session, "cutSite",        value    = "-3 18")            # Make custom PAM distance
		updateTextAreaInput(     session, "overhang",       value    = "0 5")              # Make custom PAM overhang
	})
	
	####################################
	######## EXAMPLE: Ensembl ##########
	####################################
	observeEvent(input$exampleEnsembl, {
		reset()        # Clear all the inputs
		resetOutputs() # Clear the output regions
		
		# Reset/clear the input exon table
		rValues$rhFrame   <- dfEmpty
		rValues$pamFrame  <- pamExample
		rValues$resetVal  <- TRUE
		
		# Do not display the input exon table
		geneSeqExampleFlag <<- FALSE
		
		#### Update functions ####
		updateCheckboxGroupInput(session, "casType",        selected = c("NRG", "TTTN"))       # Update the PAM selection
		updateRadioButtons(      session, "inputType",      selected = 3)                      # Update the inputType to Ensembl
		updateRadioButtons(      session, "firstExon",      selected = 1)                      # Update to use the first exon
    updateRadioButtons(      session, "exonTargetType", selected = 0)                      # Search all exons
		updateRadioButtons(      session, "customCutOpt",   selected = 1)                      # Use custom PAM input
		updateTextAreaInput(     session, "ensemblId",      value    = "ENSDART00000011520.8") # Pre-populate with the flh zebrafish transcript
		updateTextAreaInput(     session, "customPamSeq",   value    = "NG YYN")               # Make custom PAMs
		updateTextAreaInput(     session, "cutSite",        value    = "-3 18")                # Make custom PAM distance
		updateTextAreaInput(     session, "overhang",       value    = "0 5")                  # Make custom PAM overhang
	})
	
	####################################
	######## EXAMPLE: Copy/Paste #######
	####################################
	# TODO Add TALENs to this one
	observeEvent(input$exampleGeneSeq, {
		reset()        # Clear all the inputs
		resetOutputs() # Clear the output regions
		
		# Reset/clear the input exon table
		rValues$rhFrame   <- dfExample
		rValues$pamFrame  <- pamExample
		rValues$resetVal  <- TRUE
		
		# Set flag to display 'exon' information
		geneSeqExampleFlag <<- TRUE
		
		# Sample PAMs
		updateCheckboxGroupInput(session, "casType",       selected = c("NNGTGA", "TTTN"))   
		updateRadioButtons(      session, "inputType",     selected = 2)               # Change input type to copy/paste
		updateRadioButtons(      session, "pasteExonType", selected = 1)               # Change to custom exon input
		updateRadioButtons(      session, "customCutOpt",  selected = 1)               # Use custom PAM input
		updateTextAreaInput(     session, "customPamSeq",  value    = "NG YYN")        # Make custom PAMs
		updateTextAreaInput(     session, "cutSite",       value    = "-3 18")         # Make custom PAM distance
		updateTextAreaInput(     session, "overhang",      value    = "0 5")           # Make custom PAM overhang
		updateRadioButtons(      session, "talenOp",       selected = 1)               # Use TALENs
		updateNumericInput(      session, "armin",         value    = 15)              # Set minimum TALEN arm length to 15
		updateNumericInput(      session, "armax",         value    = 16)              # Set maximum TALEN arm length to 16
		updateRadioButtons(      session, "spacer",        selected = 0)               # Set spacer to 14 nts
		
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
	})
	
	####################################
	######## Reset Functions ###########
	####################################
	
	# Call the reset functions when the reset link is pressed
	observeEvent(input$reset, {
		reset()
		resetOutputs()
	})
	
	# Reset the outputs, depending on which input type is active
	resetOutputs <- function(){
		if(input$inputType == 1){
			resetGenBankOutputs()
			
		} else if(input$inputType == 2){
			resetGeneSeqOutputs()
			
		} else if(input$inputType == 3){
			resetEnsemblOutputs()
			
		}
		
		# Set the download flag to not display the download button or filtering options
		rValues$downloadF    <- FALSE
		rValues$downloadFGB  <- FALSE
		rValues$downloadFEns <- FALSE
		rValues$filtOpsGB    <- FALSE
		rValues$filtOpsGS    <- FALSE
		rValues$filtOpsE     <- FALSE
	}
	
	# Reset function
	reset <- function(){
		# Reset the geneSeqExample flag input
		geneSeqExampleFlag <<- FALSE
		
		updateCheckboxGroupInput(session, "casType",      selected = "NGG") # Clear pre-gen PAM selection
		updateTextAreaInput(     session, "customPamSeq", value    = "") # Clear custom pam input
		updateTextAreaInput(     session,	"cutSite",      value    = "") # Clear custom DSB
		updateTextAreaInput(     session, "overhang",     value    = "") # Clear custom overhangs
		updateRadioButtons(      session, "customCutOpt", selected = 0) # Update custom PAM sequence to 'no'
		
		if(input$inputType == 2){
			# Reset geneSeq input
			updateTextAreaInput(session, "geneSeq", value = "")
			
			updateRadioButtons(session, "pasteExonType", selected = 0)
			
		} else if(input$inputType == 1){
			# Reset genbank input
			updateTextAreaInput(session, "genbankId", label = "", value = "")
			
		} else if(input$inputType == 3){
			updateTextAreaInput(session, "ensemblId", label = "", value = "")
			
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
		updateRadioButtons(session, "inputType", selected = 1)
		
		# Reset talen inputs
		updateRadioButtons(session, "spacer", selected = 0)
		updateNumericInput(session, "armin",  value    = 15)
		updateNumericInput(session, "armax",  value    = 15)
		
		# Reset talen select
		updateRadioButtons(session, "talenOp", selected = 0)
		
		# Reset threshold input
		updateNumericInput(session, "threshold", value = 1.5)
		
		# Reset intronic controls
		#updateRadioButtons(session, "contextWiggleType", selected = 0)
		#updateRadioButtons(session, "gRNAWiggleType",    selected = 0)
		
		# Reset rhandsontable
		rValues$resetVal <- TRUE
		rValues$rhFrame  <- dfEmpty
		rValues$pamFrame <- pamEmpty
		
		#
		rValues$validExonFlag <- TRUE
		
		# Set the download flag to not display the download button or filtering options
		rValues$downloadF    <- FALSE
		rValues$downloadFGB  <- FALSE
		rValues$downloadFEns <- FALSE
		rValues$filtOpsGB    <- FALSE
		rValues$filtOpsGS    <- FALSE
		rValues$filtOpsE     <- FALSE
		
		# Empty the result storage
		results <<- 0
	}
	
	########################################################################################################################################################################
	################## Clear Output Areas ##################################################################################################################################
	########################################################################################################################################################################

	# Clear the GenBank unique outputs and also all output flags
	resetGenBankOutputs <- function(){
		# Clear copy/paste outputs
		output$genbankIdOutcome <- renderText({
			""
		}) 
		
		clearOutputFlags()
	}
	
	# Clear the Ensembl unique outputs and also all output flags
	resetEnsemblOutputs <- function(){
		clearOutputFlags()
		
	}
	
	# Clear the gene sequence copy/paste unique outputs and also all output flags
	resetGeneSeqOutputs <- function(){
		rValues$geneSeqError <- 0
		
		output$genbankIdOutcome <- renderText({
			""
		})
		
		clearOutputFlags()
	}
	
	# Set the output flags to clear
	clearOutputFlags <- function(){
		rValues$genbankResultsFlag <- FALSE
		rValues$geneSeqResultsFlag <- FALSE	
		rValues$ensemblResultsFlag <- FALSE
	}
	
	# Functions for updating slider input based on user actions
	observeEvent(input$armin, {
		updateSliderInput(session = session, inputId = "armax",  min = input$armin)
	})
	
	observeEvent(input$armax, {
		updateSliderInput(session = session, inputId = "armin",  max = input$armax)
	})

})

