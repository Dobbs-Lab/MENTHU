library(shiny)
library(shinyTable)
library(rhandsontable)
library(shinyIncubator)
library(stringr)
library(stringi)
library(Biostrings)

source("requiredFunctions.R")

shinyServer(function(input, output, session){
	output$exonInfo <- renderRHandsontable({
		df = data.frame(exonStart = numeric(5), exonEnd = numeric(5), stringsAsFactors = FALSE)
		rhandsontable(df) %>%
			hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE) %>%
			hot_col("exonStart", format = "0") %>%
			hot_col("exonEnd", format = "0")
	})
	
	
	####Update checkbox with select all####
	observeEvent(input$selectAll, {
		
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTN",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 #"Francisella FnCpf1: 5'-TTN-3'" = 10,
														 							 #"Francisella FnCpf1: 5'-YTN-3'" = 11,
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
														 selected = c("NGG", "NRG", "NNNRRT", "NNGRRT", "NNGTGA", "NNAGAAW", "NNNVRYAC", "NNNNGMTT"))
	})
	
	observeEvent(input$selectNone, {
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTN",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 #"Francisella FnCpf1: 5'-TTN-3'" = 10,
														 							 #"Francisella FnCpf1: 5'-YTN-3'" = 11,
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
														 selected = "NGG"
		)
	})

	#Submit Handler
	observeEvent(input$geneSeqSubmit, {
		results <- calculateMENTHU(input$casType, input$geneSeq, input$threshold)
		results <- results[order(-results$menthuScore, results$toolType, results$location),]
		print(results)
	output$geneSeqResults <- renderTable(results)
	  
	  
	})
	
	observeEvent(input$example, {
		#reset()
		updateTextAreaInput(session, "geneSeq", label = "", value = "CTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATACGCGTACCGCTAGCCAGGAAGAGTTTGTAGAAACGCAAAAAGGCCATCCGTCAGGATGGCCTTCTGCTTAGTTTGATGCCTGGCAGTTTATGGCGGGCGTCCTGCCCGCCACCCTCCGGGCCGTTGCTTCACAACGTTCAAATCCGCTCCCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTCCGACTGAGCCTTTCGTTTTATTTGATGCCTGGCAGTTCCCTACTCTCGCGTTAACGCTAGCATGGATGTTTTCCCAGTCACGACGT")
		updateCheckboxGroupInput(session, 
														 "casType", 
														 label = "",
														 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
														 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
														 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
														 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
														 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
														 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTN",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 #"Francisella FnCpf1: 5'-TTN-3'" = 10,
														 							 #"Francisella FnCpf1: 5'-YTN-3'" = 11,
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
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
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTN",
														 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
														 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
														 							 #"Francisella FnCpf1: 5'-TTN-3'" = 10,
														 							 #"Francisella FnCpf1: 5'-YTN-3'" = 11,
														 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
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