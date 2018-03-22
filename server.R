library(shiny)
library(shinyTable)
library(shinyjs)
library(rhandsontable)
library(shinyIncubator)
library(stringr)
library(stringi)
library(Biostrings)
library(rentrez)
library(rlist)
library(DT)
library(xlsx)

source("apeShiftFunctions.R")
source("genbankAccessoryFunctions.R")
source("menthu2.0AccessoryFunctions.R")
source("required2.0Functions.R")
source("targetAccessoryFunctions.R")



shinyServer(function(input, output, session){
	########################################################
	###################Global Variables#####################
	########################################################
	
	#Storage for final results so that it is accessible to 
	#download handler and submit buttons
	results <<- 0
	#Flag for displaying download button when there are
	#results available to download
	dF <<- reactiveValues(downloadF = FALSE)
	dFGB <<- reactiveValues(downloadFGB = FALSE)
	
	#Flag for displaying example table when clicking example
	#geneSeq link
	geneSeqExampleFlag <<- FALSE
	
	
	########################################################
	##################Validation Checks#####################
	########################################################
	
	####Make sure GenBank/RefSeq ID is properly formatted####
	validRefSeqGBGenbankId <- reactive({
		if(input$genbankId != ""){
			if(stringr::str_detect(input$genbankId, regex("^(NM|NR|XM|XR)_[0-9]{6}", ignore_case = TRUE))){
				#validate(
				#	need(1 == 2,
				#			 "Error: This ID matches RefSeq RNA accession format. Please submit an accession corresponding to a DNA sequence."))
				
			} else if(stringr::str_detect(input$genbankId, regex("^(AP|NP|YP|XP|WP)_[0-9]{6}", ignore_case = TRUE))){
				validate(
					need(1 == 2,
							 "Error: This ID matches RefSeq protein accession format. Please submit an accession corresponding to a DNA sequence."))
				
				
			} else if(stringr::str_detect(input$genbankId, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}", ignore_case = TRUE))){
				validate(
					need(1 == 2,
							 "Error: This ID matches a RefSeq complete genomic molecule, incomplete genomic region, contig, scaffold, or complete genome. We do not currently support any of these reference types due to issues surrounding exon identification. You can use a GenBank or RefSeq nucleotide entry corresponding to your gene of interest.")
				)
			
			} else if(stringr::str_detect(input$genbankId, regex("^[A-Z]{3}[0-9]{5}", ignore_case = TRUE))){
				validate(
					need(1 == 2,
							 "Error: This ID matches GenBank protein accession format. Please submit an accession corresponding to a DNA sequence."))
				
			} else {
				validate(
					need((((stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{2}[0-9]{6}", ignore_case = TRUE))) | 
								 	(stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{1}[0-9]{5}", ignore_case = TRUE)))) |
									(stringr::str_detect(input$genbankId, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}\\.[0-9]", ignore_case = TRUE)))) |
							 	(stringr::str_detect(input$genbankId, regex("^[a-zA-Z]{4}[0-9]{8,10}", ignore_case = TRUE))),
							 
							 "Error: The entered ID does not match a known Genbank NUCLEOTIDE or RefSeq NUCLEOTIDE ID format. Please check your submitted accession to make sure you are using a NUCLEOTIDE entry."))
			}
		} 
	})
	
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
				
				#Prevent users from blowing up the server
				need(nchar(input$geneSeq) < 5000, 
						 paste0("The DNA sequence has >5000 nucleotides. For sequences of this size,",
						 			 " please use the local version of MENTHU, which can be accessed via the 'Tools and Downloads' tab."))
			)
			
		} else if(input$inputType == 2) {
			#prevent running on empty submission
			validate(
				need(input$geneSeq != "", "")
			)
		}
	})
	
	#Valid threshold - make sure it's not negative
	validThreshold <- reactive({
		if(!is.null(input$threshold)){
			validate(
				need(input$threshold >= 0,
						 "Error: Threshold value must be non-negative.")
			)
		}	
	})
	
	#Valid PAM
	validPAM <- reactive({
		validate(
			need((length(input$casType) > 0) | (input$talenOp == 1) | input$customCutOpt == 1,
					 "Error: No nuclease selected. Please select a Cas type, choose to use a custom PAM scheme, and/or the TALEN input option.")
		)
	})
	
	
	#Valid TALEN inputs
	validTalen <- reactive({
		validate(
			need(input$armin <= input$armax,
					 "Error: Maximum TALEN arm length must be greater than or equal to minimum TALEN arm length."),
			
			need(input$spamin <= input$spamax,
					 "Error: Maximum spacer length must be greater than or equal to minimum spacer length."),
			
			need((input$spamin >= 14),
					 "Error: Minimum spacer length is 14 nucleotides."),
			
			need((input$spamax <= 16),
					 "Error: Maximum spacer length is 16 nucleotides."),
			
			need((input$armin >= 15),
					 "Error: Minimum TALEN arm length is 15 nucleotides."),
			
			need((input$armax <= 18),
					 "Error: Maximum TALEN arm length is 18 nucleotides.")
			
		)
	})
	
	#Valid custom PAM inputs
	validCustomPam <- reactive({
		#If the custom PAM is selected, make sure it only has IUPAC nucleotides or separating characters
		if(input$customPamSeq != ""){
			if(stringr::str_detect(input$customPamSeq, "[^ACGTRYSWKMBDHVNacgtryswkmbdhvn\\s,]")){
				validate(
					need(1 == 2,
							 "Error: Non-allowed characters detected. Only standard nucleotide (A, C, G, T), IUPAC extended nucleotide symbols (R, Y, S, W, K, M, B, D, H, V, N), and separating characters (spaces and commas) allowed.")
				)
			} else {
				disallowed <- "^[N]+[ACGTRYSWKMBDHV]{1}$|^[N]+[ACGTRYSWKMBDHV]{1}[N]+$|^[ACGTRYSWKMBDHVN]{1}$|^[N]+$"
				disallowedCheck <- grepl(disallowed, pamStitch("", input$customPamSeq), ignore.case = TRUE, perl = TRUE)
				disallowedExist <- is.element(TRUE, unlist(disallowedCheck))
				if(disallowedExist){
					validate(
						need(1 == 2,
								 "Error: Due to computational limitations, we do not accept custom PAMs consisting solely of 'N', single nucleotide PAMs, or PAMs consisting of solely of 'N's and a single nucleotide.")
					)
				}
			}
			
		}
	})
	
	#Valid custom cut sites
	validCustomCutSites <- reactive({
		if(input$cutSite != ""){
			validate(
				need(!stringr::str_detect(input$cutSite, "[^0-9\\s,-]"),
						 "Error: Only numbers, dashes, commas, and spaces allowed.")
			)
		}
	})
	
	#Check that the number of custom PAMs matches the number of custom cutSites
	validMatchCustomInputLength <- reactive({
		if((input$customPamSeq != "") && 
			 (input$customCutOpt != "")){
			validate(
				need(length(as.character(distStitch("", input$cutSite))) == length(pamStitch("", input$customPamSeq)),
				"Error: The number of custom PAM sequences does not match the number of custom DSB sites. Please make sure that each PAM has one distance to its cut site specified.")
			)
		}
	})
	
	########################################################
	##############PRINT VALIDATION RESULTS##################
	########################################################
	output$validgenbankid <- renderText({
		#validGenBankId()
		validRefSeqGBGenbankId()
	})
	
	output$genbankidexists <- renderText({
		
	})
	
	output$validgeneseq <- renderText({
		validGeneSeq()
	})
	
	output$validthreshold <- renderText({
		validThreshold()
	})
	
	output$validpam <- renderText({
		validPAM()
	})
	
	output$validtalen <- renderText({
		validTalen()
	})
	
	output$validcustompam <- renderText({
		validCustomPam()
	})
	
	output$validcustomcutsites <- renderText({
		validCustomCutSites()
	})
	
	output$validmatchcustominputlength <- renderText({
			validMatchCustomInputLength()
	})
	
	#output$validexoninfo <- renderText({
	#	validExonInfo()
	#})
	
	########################################################
	#####################UI Rendering#######################
	########################################################
	
	#Render exon table
	output$exonInfo <- renderRHandsontable({
		exonOutDF = renderTable()
		if(!is.null(exonOutDF)){
			rhandsontable(exonOutDF) %>%
				hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>% 
				hot_col("exonStart", format = "0") %>%
				hot_col("exonEnd",   format = "0")
		}
	})
	
	#render function for RHandsonTable
	renderTable <- reactive({
		if(!is.null(input$exonInfo)){
			exonDF = hot_to_r(input$exonInfo)
		} else {
			if(geneSeqExampleFlag){
				exonDF = data.frame(exonStart        = c(1, 101, 201, 301, 401),
														exonEnd          = c(100, 200, 300, 400, 500),
														stringsAsFactors = FALSE)
			} else {
				exonDF = data.frame(exonStart        = numeric(5), 
														exonEnd          = numeric(5), 
														stringsAsFactors = FALSE)
			}
		}
	})
	
	#Download button
	output$downOut <- renderUI({
		if(dF$downloadF){
			downloadButton("downRes", "Download Results")
		} else {
			""
		}
	})
	
	output$downOutGenbank <- renderUI({
		if(dFGB$downloadFGB){
			downloadButton("downRes", "Download Results")
		} else {
			""
		}
	})
	
	#Download handler
	output$downRes <- downloadHandler(
		filename = function(){
			#Name file in the form "YYYY-MM-DD_HH-MM-SS_targets.csv
			paste(gsub("CDT", "", gsub(" ", "_", Sys.time())), "_targets.csv")},
		
		content = function(file){
			write.csv(results, file, row.names = FALSE)
		}
		
	)
	
	
	########################################################
	#################Submission Handling####################
	########################################################
	
	observeEvent(input$genbankSubmit,{
		talFlag <- 2
		#Run checks for okay PAM/TALEN input
		if(input$talenOp == 1){
			if(is.null(validTalen())){
				talFlag <- 1 #Input TALEN is good
			} else {
				talFlag <- 2 #Problems with input TALEN
			}
		} else {
			talFlag <- 0 #TALENs not selected
		}
		
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
		
		#Prevent from running without okay inputs
		if(is.null(validRefSeqGBGenbankId()) && #Check Genbank ID is okay
			 is.null(validThreshold()) &&         #Check threshold is okay
			 is.null(validPAM()) &&               #Check that one of the input options is selected
			 (cusPamFlag != 2) &&                 #Check if custom PAMs are used, and if they are okay
			 (talFlag != 2) &&
			 (lenMatch != 2)){                     #Check if TALENs are used, and if they are okay
			
			progress <- shiny::Progress$new()
			#Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			#Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			talenList <- ""
			
			if(input$talenOp == 1){
				talenList <- list(input$armin, input$armax, input$spamin, input$spamax)
			}
			
			progress$set(detail = "Retrieving GenBank entry...", value = 0.1)
			
			#Try to pull genbank entry associated with accession
			#Get the GenBank sequence with exon/intron information
			gba <- input$genbankId
			
			endOfTry <<- FALSE
			gbFlag   <<- FALSE
			gbhFlag  <<- FALSE
			
			#Try to retrieve the Genbank file 
			if(endOfTry == FALSE){
				tryCatch({
					info <- suppressWarnings(wonkyGenBankHandler(gba))
					endOfTry <<- TRUE
					gbhFlag  <<- TRUE  #Flag to indicate wonkyGenBankHandler was required 
					gbFlag   <<- FALSE #Flag to indicate readGenBank failed #deprecated; to remove
				}, error = function(err){
				}
				)
			}
			
			#Figure out which exons to target
			if(endOfTry){
				if(input$exonTargetType == 0){
					exonStuff <- 1
				} else	if(input$exonTargetType == 1){
					exonStuff <- input$exonBegPercentage/100
				} else if(input$exonTargetType == 2){
					exonStuff <- input$exonEndPercentage/100
				} else if(input$exonTargetType == 3){
					exonStuff <- input$exonTargetList	
				}
				
				#Handle cut distances and PAMs for input to calculateMENTHUGeneSeqGenBank
				if(input$customCutOpt == 1){ #If customs pams in use
					if(input$casType != ""){ #If pre-made PAMs in use
						pams <- pamStitch(input$casType, input$customPamSeq)
						cutDistances <- distStitch(input$casType, input$cutSite)
					} else { #If no pre-made PAMs
						pams <- pamStitch("", input$customPamSeq)
						cutDistances <- distStitch("", input$cutSite)
					}
					
				} else { #If no custom pams
					pams <- input$casType
					cutDistances <- rep(-3, length(input$casType))
				}
				
				#Calculate the MENTHU score
				results <<- calculateMENTHUGeneSeqGenBank(pams, cutDistances, wiggle = TRUE, wigRoom = 39, talenList, gbFlag, gbhFlag, 
																									info, input$threshold, input$firstExon, input$exonTargetType, exonStuff, progress, input$scoreScheme)
				#Order the result table from largest menthuScore to smallest
				results <<- results[order(-results$MENTHU_Score),]

				#Set the download button flag to true to render download button visible
				dFGB$downloadFGB <<- TRUE
				
				#Output the results in a data table
				output$genbankResults <- DT::renderDataTable(results, 
																										 options = list(scrollX = TRUE), 
																										 rownames = FALSE)
				
			} else {
				#Output error if no genbank file is found
				output$genbankIdOutcome <- renderText(paste0("Error: Accession '", input$genbankId, "' was not found in database."))
			}
		}
	})
	
	#Copy/paste Submit Handler
	observeEvent(input$geneSeqSubmit, {
		
		#Run checks for okay PAM/TALEN input
		talFlag <- 2
		if(input$talenOp == 1){
			if(is.null(validTalen())){
				talFlag <- 1 #Input TALEN is good
			} else {
				talFlag <- 2 #Problems with input TALEN
			}
		} else {
			talFlag <- 0 #TALENs not selected
		}
		
		#Run checks for okay custom PAM input
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
		
		#Run checks for okay length match b/w PAM and cut site
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
		
		#Prevent from running without okay inputs
		if(is.null(validGeneSeq()) && #Check Genbank ID is okay
			 is.null(validThreshold()) &&         #Check threshold is okay
			 is.null(validPAM()) &&               #Check that one of the input options is selected
			 (cusPamFlag != 2) &&                 #Check if custom PAMs are used, and if they are okay
			 (talFlag != 2) &&
			 (lenMatch != 2)){                     #Check if TALENs are used, and if they are okay
			
			#Make a progress object
			progress <- shiny::Progress$new()
			#Make sure the progress option closes regardless of calculation outcome
			on.exit(progress$close())
			
			#If there is not exon input, just set to 0, otherwise use exon input
			if(input$pasteExonType == 0){
				exonIn <- 0
			} else {
				if(is.null(input$exonInfo)){
					exonIn <- 0
				} else {
					exonIn <- exonHandler(input$exonInfo)
				}
				
			}
			
			#Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			if(input$talenOp == 0){
				talArmin <- ""
				talArmax <- ""
				talSpamin <- ""
				talSpamax <- ""
			} else {
				talArmin <- input$armin
				talArmax <- input$armax
				talSpamin <- input$spamin
				talSpamax <- input$spamax
			}
			
			#Format custom PAMs and cut distances for submission to calculateMENTHUGeneSeq
			if(input$customCutOpt == 1){ #If customs pams in use
				if(input$casType != ""){ #If pre-made PAMs in use
					pams <- pamStitch(input$casType, input$customPamSeq)
					cutDistances <- distStitch(input$casType, input$cutSite)
				} else { #If no pre-made PAMs
					pams <- pamStitch("", input$customPamSeq)
					cutDistances <- distStitch("", input$cutSite)
				}
			} else { #If no custom pams
				pams <- input$casType
				cutDistances <- rep(-3, length(input$casType))
			}
			
			#Calculate the MENTHU score
			results <<- calculateMENTHUGeneSeq(pams, cutDistances, wiggle = TRUE, wigRoom = 39, input$geneSeq, input$threshold, 
																				 exonIn, progress, talArmin, talArmax, talSpamin, talSpamax, input$scoreScheme)
			#Order the result table from largest menthuScore to smallest
			results <<- results[order(-results$MENTHU_Score),]

			#Set the download button flag to true to render download button visible
			dF$downloadF <<- TRUE
			#Output the results in a data table
			output$geneSeqResults <- DT::renderDataTable(results, 
																									 options = list(scrollX = TRUE), 
																									 rownames = FALSE)
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
	
	observeEvent(input$exampleGenBank, {
		reset()
		geneSeqExampleFlag <<- FALSE
		updateRadioButtons(session,
											 "inputType", 
											 selected = 1)
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NRG", "NNNRRT"))
		updateTextAreaInput(session,
												"genbankId",
												value = "AY214391.1")
		updateRadioButtons(session,
											 "firstExon",
											 selected = 0)
		updateRadioButtons(session,
											 "exonTargetType",
											 selected = 0)
		if(input$scoreScheme == 1){
			updateNumericInput(session,
												 "threshold",
												 value = 40)
		} else {
			updateNumericInput(session,
												 "threshold",
												 value = 1.5)
		}
		
		#Use custom PAM input
		updateRadioButtons(session,
											 "customCutOpt",
											 selected = 1)
		
		#Make custom PAMs
		updateTextAreaInput(session,
												"customPamSeq",
												value = "NNNGCT")
		
		#Make custom PAM distance
		updateTextAreaInput(session,
												"cutSite",
												value = "-4")
		
		
	})
	
	#Example of gene sequence copy/paste input
	observeEvent(input$exampleGeneSeq, {
		reset()
		
		#Set flag to display 'exon' information
		geneSeqExampleFlag <<- TRUE
		
		#Change input type to copy/paste
		updateRadioButtons(session,
											 "inputType",
											 selected = 2)
		
		#Change to custom exon input
		updateRadioButtons(session,
											 "pasteExonType",
											 selected = 1)
		
		#Sample copy/paste gene input
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
		
		#Sample PAMs
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
		
		#Use custom PAM input
		updateRadioButtons(session,
											 "customCutOpt",
											 selected = 1)
		
		#Make custom PAMs
		updateTextAreaInput(session,
												"customPamSeq",
												value = "NNNGAG NNNGNG")
		
		#Make custom PAM distance
		updateTextAreaInput(session,
												"cutSite",
												value = "-4 -3")
		
		
		#Change threshold to 30 to show output
		updateNumericInput(session, "threshold", value = 30)
		
		#
		

	})
	
	#When hitting the reset link button, call the reset function
	observeEvent(input$reset, {
		reset()
		resetOutputs()
	})
	
	#Reset function
	reset <- function(){
		#Reset the geneSeqExample flag input
		geneSeqExampleFlag <<- FALSE
		
		#Update checkboxes to only have SpCas9 NGG selected
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = "NGG"
		)
		
		#Clear custom pam input
		updateTextAreaInput(session,
												"customPamSeq",
												value = ""
		)
		
		#Clear custom DSB
		updateTextAreaInput(session,
												"cutSite",
												value = "")
		
		#Update custom PAM sequence to 'no'
		updateRadioButtons(session,
											 "customCutOpt",
											 selected = 0)

		#Reset geneSeq input
		updateTextAreaInput(session, 
												"geneSeq", 
												label = "", 
												value = "")
		
		#Reset genbank input
		updateTextAreaInput(session,
												"genbankId",
												label = "",
												value = "")
		
		#Reset sequence input type
		updateRadioButtons(session, 
											 "inputType", 
											 selected = 1
		)
		
		
		updateRadioButtons(session,
											 "pasteExonType",
											 selected = 0)
		
		#Reset talen inputs
		updateNumericInput(session, "spamin", value = 14)
		updateNumericInput(session, "spamax", value = 16)
		updateNumericInput(session, "armin",  value = 15)
		updateNumericInput(session, "armax",  value = 15)
		
		#Reset talen select
		updateRadioButtons(session, "talenOp", selected = 0)
		
		#Reset threshold input
		updateNumericInput(session, "threshold", value = 40)
		
		#Empty the result storage
		results <<- 0
		#Set the download flag to not display the download button
		dF$downloadF <- FALSE
		dFGB$downloadFGB <- FALSE
		
	}
	
	resetOutputs <- function(){
		#Clear outputs
		output$geneSeqResults <- renderText({
			""
		})
		
		output$genbankIdOutcome <- renderText({
			""
		})

	}
	
	#Functions for updating slider input based on user actions
	observeEvent(input$armin, {
		updateSliderInput(session = session, inputId = "armax", min = input$armin)
	})
	
	observeEvent(input$armax, {
		updateSliderInput(session = session, inputId = "armin", max = input$armax)
	})
	
	observeEvent(input$spamin, {
		updateSliderInput(session = session, inputId = "spamax", min = input$spamin)
	})
	
	observeEvent(input$spamax, {
		updateSliderInput(session = session, inputId = "spamin", max = input$spamax)
	})
	
})