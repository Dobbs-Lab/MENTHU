library(shiny)
library(shinyTable)
library(rhandsontable)
library(shinyIncubator)
library(stringr)
library(stringi)
library(Biostrings)
library(genbankr)
library(DT)
library(xlsx)

source("requiredFunctions.R")
source("apeShiftFunctions.R")
source("targetAccessoryFunctions.R")
source("menthuAccessoryFunctions.R")
source("genbankAccessoryFunctions.R")

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
	validGenBankId <- reactive({
		#If input type is genbank/refseq, and the string is not empty
		if((input$inputType == 1) && (input$genbankId != "")){
			
			#If the input starts off with two characters followed by an underscore
			if(str_detect(toupper(input$genbankId), "[A-Za-z0-9]{2}\\_")){
				
				#If the input has more than one underscore
				if(str_detect(toupper(input$genbankId), "[\\_].*[\\_]")){
					validate(
						need(!str_detect(toupper(input$genbankId), "[\\_].*[\\_]"),
								 "Error: This accession contains two underscores. RefSeq accession IDs contain only a single underscore.")
					)
					
				} else {
					validate(
						#If the sequence contains anything other than alphanumerics, underscores, or periods
						need(!str_detect(toupper(input$genbankId), "[^\\w\\.]"),
								 paste0("Error: The accession contains restricted characters. ",
								 			 "Accession IDs may contain only letters, numbers, periods,",
								 			 " and a single underscore ('_') after the first two characters, in the case of RefSeq IDs.")
						))
				}
				
			} else {
				#If a genbank accession, check that it follows genbank format
				validate(
					need(!str_detect(toupper(input$genbankId), "[^A-Za-z0-9\\.]"),
							 paste0("Error: The accession contains restricted characters.",
							 			 " Accession IDs may contain only letters, numbers, periods, ",
							 			 "and a single underscore ('_') after the first two characters, in the case of RefSeq IDs.")
					)
				)
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
				need(!str_detect(input$geneSeq, "[>]"), 
						 "Error: Input DNA sequence appears to be in FASTA format. Please paste your sequence without the fasta header."),
				
				#Check for not DNA input
				need(!str_detect(input$geneSeq, "[^ACGTacgt0-9\\s\\n]"), 
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
	
	#Exon input validation
	#validExonInfo <- reactive({
	#	if(!is.null(input$exonInfo)){
	#		values[["old"]] <- isolate(values[[""]])
	#	}
	#})
	
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
		#print(input$casType) #For testing purposes
		validate(
			need((length(input$casType) > 0) | (input$talenOp == 1),
					 "Error: No PAMs selected.")
		)
	})
	
	
	########################################################
	#####################UI Rendering#######################
	########################################################
	
	#Render exon table
	output$exonInfo <- renderRHandsontable({
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
		
		rhandsontable(exonDF) %>%
			hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>% 
			hot_col("exonStart", format = "0") %>%
			hot_col("exonEnd",   format = "0")
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
			write.csv(results, file)
		}
		
	)
	
	########################################################
	##############PRINT VALIDATION RESULTS##################
	########################################################
	output$validgenbankid <- renderText({
		validGenBankId()
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
	
	#output$validexoninfo <- renderText({
	#	validExonInfo()
	#})
	
	###Exon Input Table Handler####
	#observe({
	#	if(!is.null(input$exonInfo)){
	#		values[["prev"]] <- isolate(values[["exonDF"]])
	#		exonDF = hot_to_r(input$exonInfo)
	#	} else {
	#		if(is.null(values[["exonDF"]])){
	#			exonDF <- exonDF
	#		} else {
	#			exonDF <- values[["exonDF"]]
	#		}
	#	}
	#	values[["exonDF"]] <- exonDF
	#})
	
	observeEvent(input$genbankSubmit,{
		progress <- shiny::Progress$new()
		#Make sure the progress option closes regardless of calculation outcome
		on.exit(progress$close())
		
		
		#Prevent from running if missing inputs or have invalid inputs
		#if(is.null(validThreshold()) && 
		#	 (is.null(validPAM()))){
			
			#Get exon info
			
			#Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			#exonPer <- 0 #Dummy variable for calculateMENTHUGeneSeq function call
			
			talenList <- ""
			
			if(input$talenOp == 1){
				talenList <- list(input$armin, input$armax, input$spamin, input$spamax)
			}
			
			progress$set(detail = "Retrieving GenBank entry...", value = 0.1)
			
			#Try to pull genbank entry associated with accession
			#Get the GenBank sequence with exon/intron information
			#Make a genbank acession object
			gba <- genbankr:::GBAccession(input$genbankId)
			
			endOfTry <<- FALSE
			gbFlag   <<- FALSE
			gbhFlag  <<- FALSE
			
			#Experimental try-catch to deal with exon/cds files
			tryCatch({
				info <- genbankr:::readGenBank(gba, partial = TRUE, verbose = TRUE)
				
				endOfTry <<- TRUE
				gbFlag   <<- TRUE #Flag to indicate readGenBank succeeded
				
			}, error = function(err){
				endOfTry <<- FALSE
				
			}
			)
			
			print(endOfTry)
			
			if(endOfTry == FALSE){
				tryCatch({
					#print("We got to this point") #Bughunting
					info <- suppressWarnings(wonkyGenBankHandler(gba))
					
					endOfTry <<- TRUE
					gbhFlag  <<- TRUE  #Flag to indicate wonkyGenBankHandler was required
					gbFlag   <<- FALSE #Flag to indicate readGenBank failed
				}, error = function(err){
					output$validgenbankid <- renderText({
						validGenBankId()
					})
				}
				)
			}
			
			print(endOfTry)
			
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
				print(input$inputType)
				#Calculate the MENTHU score
				results <<- calculateMENTHUGeneSeqGenBank(input$casType, wiggle = TRUE, wigRoom = 39, talenList, gbFlag, gbhFlag, info, input$threshold, input$firstExon, input$exonTargetType, exonStuff, progress)
				#Order the result table from largest menthuScore to smallest
				results <<- results[order(-results$MENTHU_Score),]
				#print(results)
				#Set the download button flag to true to render download button visible
				dFGB$downloadFGB <<- TRUE
				#Output the results in a data table
				output$genbankResults <- DT::renderDataTable(results, 
																										 options = list(scrollX = TRUE), 
																										 rownames = FALSE)
				
			}

		#}
	})
	
	#Submit Handler
	observeEvent(input$geneSeqSubmit, {
		#Make a progress object
		progress <- shiny::Progress$new()
		#Make sure the progress option closes regardless of calculation outcome
		on.exit(progress$close())
		
		
		#Prevent from running if missing inputs or have invalid inputs
		if(is.null(validGeneSeq()) && 
			 (is.null(validThreshold())) && 
			 (is.null(validPAM()))){
			
			#If there is not exon input, just set to 0, otherwise use exon input
			if(input$pasteExonType == 0){
				exonIn <- 0
			} else {
				exonIn <- exonHandler(input$exonInfo)
				if(nrow(exonIn) == 0){
					exonIn <- 0
				}
			}
			
			print(exonIn)
			
			#Set the progress message to display at beginning of calculations
			progress$set(message = "Progress:", value = 0)
			
			#results <- calculateMENTHUGeneSeq(input$casType, input$geneSeq, input$threshold, exonIn)
			#	exonList <- 0
			#} else if(input$exonTargetType == 1){
			#	exonList <- -(input$exonBegPercentage)
			#} else if(input$exonTargetType == 2){
			#	exonList <- input$exonEndPercentage
			#} else if(input$exonTargetType == 3){
			#	exonList <- input$exonTargetList
			#}
			#exonPer <- 0 #Dummy variable for calculateMENTHUGeneSeq function call
			
			#talenList <- ""
			
			#if(input$talenOp == 1){
			#	talenList <- list(input$armin, input$armax, input$spamin, input$spamax)
			#}
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
			

			#Calculate the MENTHU score
			results <<- calculateMENTHUGeneSeq(input$casType, wiggle = TRUE, wigRoom = 39, input$geneSeq, input$threshold, exonIn, progress, talArmin, talArmax, talSpamin, talSpamax)
			#Order the result table from largest menthuScore to smallest
			results <<- results[order(-results$MENTHU_Score),]
			#print(results)
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
		updateRadioButtons(session,
											 inputId = "talenOp",
											 selected = 1)
	})
	
	observeEvent(input$exampleGenBank, {
		reset()
		geneSeqExampleFlag <<- FALSE
		updateRadioButtons(session,
											 "inputType", 
											 selected = 1)
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NGG"))
		updateTextAreaInput(session,
												"genbankId",
												value = "U49845")
		updateRadioButtons(session,
											 "firstExon",
											 selected = 1)
		updateRadioButtons(session,
											 "exonTargetType",
											 selected = 3)
		updateTextAreaInput(session,
												"exonTargetList",
												value = "1,3")
		updateNumericInput(session,
											 "threshold",
											 value = 30)
	})
	
	#Example of gene sequence copy/paste input
	observeEvent(input$exampleGeneSeq, {
		reset()
		geneSeqExampleFlag <<- TRUE
		
		updateRadioButtons(session,
											 "inputType",
											 selected = 2)
		
		updateRadioButtons(session,
											 "pasteExonType",
											 selected = 1)
		
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
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = c("NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
		updateNumericInput(session, "threshold", value = 30)
		
	})
	
	#When hitting the reset link button, call the reset function
	observeEvent(input$reset, {
		reset()
		print(dF$downloadF)
		resetOutputs()
	})
	
	#Reset function
	reset <- function(){
		#Update checkboxes to only have SpCas9 NGG selected
		updateCheckboxGroupInput(session, 
														 "casType", 
														 selected = "NGG"
		)
		
		#Reset geneSeq input
		updateTextAreaInput(session, "geneSeq", label = "", value = "")
		
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
		
		
		
		
		
		#Reset threshold input
		updateNumericInput(session, "threshold", value = 40)
		
		#Empty the result storage
		results <<- 0
		#Set the download flag to not display the download button
		dF$downloadF <- FALSE
		dFGB$downloadFGB <- FALSE
		#Reset the geneSeqExample flag input
		geneSeqExampleFlag <<- FALSE
	}
	
	resetOutputs <- function(){
		#Clear outputs
		output$geneSeqResults <- renderText({
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