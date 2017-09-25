library(shiny)
library(shinyTable)
library(rhandsontable)
library(shinyIncubator)
library(stringr)
library(stringi)
library(Biostrings)
library(DT)

source("requiredFunctions.R")

shinyServer(function(input, output, session){
	####Validation Checks####
	
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
								 "Error: The accession contains restricted characters. Accession IDs may contain only letters, numbers, periods, and a single underscore ('_') after the first two characters, in the case of RefSeq IDs."
						))
				}
				
			} else {
				#If a genbank accession
				validate(
					need(!str_detect(toupper(input$genbankId), "[^A-Za-z0-9\\.]"),
							 "Error: The accession contains restricted characters. Accession IDs may contain only letters, numbers, periods, and a single underscore ('_') after the first two characters, in the case of RefSeq IDs."
					)
				)
			}
			
		} 
	})
	
	####Validate copy/paste sequence input####
	validGeneSeq <- reactive({
		#If the input type is copy/paste gene seq and text has been entered
		if((input$geneSeq != "") && (input$inputType == 2)){
			validate(
				#Check for fasta format
				need(!str_detect(input$geneSeq, "[>]"), 
						 "Error: Input DNA sequence appears to be in FASTA format. Please paste your sequence without the fasta header."),
				#Check for not DNA input
				need(!str_detect(input$geneSeq, "[^ACGTacgt0-9\\s\\n]"), 
						 paste0("Error: Input DNA sequence contains non-standard nucleotides. ", 
						 			 "Allowed nucleotides are A, C, G, and T.")),
				need(nchar(input$geneSeq) < 5000, "The DNA sequence has >5000 nucleotides. For sequences of this size, please use the local version of MENTHU, which can be accessed via the 'Tools and Downloads' tab.")
			)
		} else if(input$inputType == 2) {
			#prevent running on empty submission
			validate(
				need(input$geneSeq != "", "")
			)
		}
	})
	
	#Exon input validation
	validExonInfo <- reactive({
		if(!is.null(input$exonInfo)){
			values[["old"]] <- isolate(values[[""]])
		}
	})
	
	#Valid threshold
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
		print(input$casType)
		validate(
			need(length(input$casType) > 0,
					 "Error: No PAMs selected.")
		)
	})
	
	#Render exon table
	output$exonInfo <- renderRHandsontable({
		if(!is.null(input$exonInfo)){
			exonDF = hot_to_r(input$exonInfo)
		} else {
			exonDF = data.frame(exonStart = numeric(5), exonEnd = numeric(5), stringsAsFactors = FALSE)
		}
		
			rhandsontable(exonDF) %>%
				hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>%
				hot_col("exonStart", format = "0") %>%
				hot_col("exonEnd", format = "0")
		
	})
	
	
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
	
	####Update checkbox with select all####
	observeEvent(input$selectAll, {
		
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "1. Select the PAM sequence(s) you wish to target:",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTTN",
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 #"Francisella FnCpf1: 5'-TTN-3'" = "TTN",
														 #"Francisella FnCpf1: 5'-YTN-3'" = "YTN"),
														 selected = c("NGG", "NRG", "NNNRRT", "NNGRRT", "NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
	})
	
	observeEvent(input$selectNone, {
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "1. Select the PAM sequence(s) you wish to target:",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTTN",
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 #"Francisella FnCpf1: 5'-TTN-3'" = "TTN",
														 #"Francisella FnCpf1: 5'-YTN-3'" = "YTN"),
														 selected = "NGG"
		)
	})

	#Submit Handler
	observeEvent(input$geneSeqSubmit, {
		#Make a progress object
		progress <- shiny::Progress$new()
		#Make sure this closes
		on.exit(progress$close())
		
		
		#Prevent from running if missing inputs or have invalid inputs
		if(is.null(validGeneSeq()) && (is.null(validThreshold())) && (is.null(validPAM()))){

			if(input$pasteExonType == 0){
				exons <- 0
			} else {
				exons <- exonHandler(input$exonInfo)
				if(nrow(exons) == 0){
					exons <- 0
				}
			}
			
				print(exons)
				progress$set(message = "Beginning calculation...", value = 0)
				
				#results <- calculateMENTHUGeneSeq(input$casType, input$geneSeq, input$threshold, exons)
				results <- calculateMENTHUGeneSeq(input$casType, input$geneSeq, input$threshold, exons, progress)
				results <- results[order(-results$menthuScore),]
				#print(results)
				output$geneSeqResults <- DT::renderDataTable(results, options = list(scrollX = TRUE))

			
		}
	})

	
	observeEvent(input$exampleGeneSeq, {
		#reset()
		updateRadioButtons(session,
											 "inputType",
											 label = "2. What type of sequence input do you want to use?",
											 choices = list("GenBank Gene ID"  = 1,
											 							 "Copy/paste gene sequence" = 2),
											 selected = 2)
		updateTextAreaInput(session, "geneSeq", label = "", value = "CTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATACGCGTACCGCTAGCCAGGAAGAGTTTGTAGAAACGCAAAAAGGCCATCCGTCAGGATGGCCTTCTGCTTAGTTTGATGCCTGGCAGTTTATGGCGGGCGTCCTGCCCGCCACCCTCCGGGCCGTTGCTTCACAACGTTCAAATCCGCTCCCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTCCGACTGAGCCTTTCGTTTTATTTGATGCCTGGCAGTTCCCTACTCTCGCGTTAACGCTAGCATGGATGTTTTCCCAGTCACGACGT")
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "1. Select the PAM sequence(s) you wish to target:",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTTN",
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 #"Francisella FnCpf1: 5'-TTN-3'" = "TTN",
														 #"Francisella FnCpf1: 5'-YTN-3'" = "YTN"),
														 selected = c("NGG", "NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
		updateNumericInput(session, "threshold", value = 30)
		
	})
	
	observeEvent(input$reset, {
		reset()
	})
	
	reset <- function(){
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTTN",
														 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 #"Francisella FnCpf1: 5'-TTN-3'" = "TTN",
														 #"Francisella FnCpf1: 5'-YTN-3'" = "YTN"),
														 selected = "NGG"
		)
		updateRadioButtons(session, "inputType", label = "", choices = list(#"GenBank Gene ID"  = 1,
			"Copy/paste sequence" = 2
		))
		updateTextAreaInput(session, "geneSeq", label = "", value = "")
		updateNumericInput(session, "threshold", label = "", value = 40)
		resetOutputs()

	}
	
	resetOutputs <- function(){
		output$geneSeqResults <- renderTable({
			""
		})
	}
})