library(shiny)
library(shinyjs)
library(rhandsontable)

shinyUI(function(request){
	
	####Creates the navbar set-up####
	navbarPage(id = 'mainPage',
						 windowTitle = "MENTHU",
						 
						 #Stylesheet
						 theme = "ogtheme.css", 
						 
						 #Page title box
						 tags$div(""),
						 
						 ########ABOUT TAB#################################################
						 tabPanel(#id = 'about',
						 				 tags$div("MENTHU v2.1.0"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 				 											 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 				 											 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 				 											 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 				 											 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://github.com/Dobbs-Lab/MENTHU", 
						 				 											 target = "_blank", tags$img(src = "GitHub_Logo.png",                   width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://hub.docker.com/r/cmmann/menthu", 
						 				 											 target = "_blank", tags$img(src = "Docker_Logo.png",                   width = "100%"))))
						 				 	
						 				 )),
						 				 
						 				 #Text area in center of page
						 				 column(9, wellPanel(
						 				 	
						 				 	#Display about page
						 				 	includeHTML("www/menthuAbout.html")
						 				 ))
						 				 
						 ),
						 
						 ##########INSTRUCTIONS############################################
						 tabPanel(tags$div("Instructions and FAQs"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 				 											 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 				 											 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 				 											 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 				 											 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 				 )),
						 				 
						 				 #Text area in center of page
						 				 column(9, wellPanel(
						 				 	
						 				 	#Display page instructions
						 				 	includeHTML("www/menthuInstructions.html")
						 				 ))
						 				 
						 ),
						 
						 ##########Calculate MENTHU TAB###################################
						 tabPanel(id = "single",
						 				 tags$div("Submit Job"),
						 				 titlePanel(""),
						 				 
						 				 ##Sidebar############################################################
						 				 #Adds a sidebar for users to pre-populate fields with an example, and reset the form
						 				 column(2, wellPanel(
						 				 	class = "examplePanel",
						 				 	
						 				 	# Attempts to get the sidebar floating have failed thus far
						 				 	#style = "position:fixed;width:inherit;",
						 				 	
						 				 	p(tags$b(tags$u("Example Inputs"))),
						 				 	
						 				 	#GenBank Example
						 				 	actionLink("exampleGenBank",
						 				 						 label = "[GenBank Gene ID Example]"),
						 				 	
						 				 	tags$br(),
						 				 	tags$br(),
						 				 	
						 				 	# Ensembl input example; input$exampleEnsembl
						 				 	actionLink("exampleEnsembl",
						 				 						 label = "[Ensembl ID Example]"),
						 				 	
						 				 	tags$br(),
						 				 	tags$br(),
						 				 	
						 				 	#Cut/Paste cDNA example; input$example
						 				 	actionLink("exampleGeneSeq", 
						 				 						 label = "[Pasted Sequence Example]"),
						 				 	
						 				 	tags$br(),
						 				 	tags$br(),
						 				 	
						 				 	#Reset Button; input$reset
						 				 	actionLink("reset", 
						 				 						 label = "Reset Form")
						 				 	
						 				 )),
						 				 
						 				 ####Main Bar#########################################################
						 				 # Main panel for entering information and submitting job
						 				 column(9, 
						 				 			 
						 				 			 ####Choose PAM sequence#############################################
						 				 			 wellPanel(
						 				 			 	checkboxGroupInput("casType",
						 				 			 										 label = "1. Select the PAM sequence(s) you wish to target:",
						 				 			 										 choices = list("S. pyogenes SpCas9: 5'-NGG-3'"         = "NGG",
						 				 			 										 							 "S. pyogenes SpCas9: 5'-NRG-3'"          = "NRG",
						 				 			 										 							 "S. aureus SaCas9: 5'-NNNRRT-3'"         = "NNNRRT",
						 				 			 										 							 "S. aureus SaCas9: 5'-NNGRRT-3'"         = "NNGRRT",
						 				 			 										 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'"   = "NNGTGA",
						 				 			 										 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'"  = "NNAGAAW",
						 				 			 										 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'"       = "NNNVRYAC",
						 				 			 										 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT",
						 				 			 										 							 "Acidaminococcus AsCas12a (AsCpf1)/Lachnospiraceae LbCas12a (LbCpf1): 5'-TTTN-3'" = "TTTN",
						 				 			 										 							 "Acidaminococcus AsCas12a (AsCpf1)/Lachnospiraceae LbCas12a (LbCpf1): 5'-TTTV-3'" = "TTTV",
						 				 			 										 							 "Francisella FnCas12a (FnCpf1): 5'-TTN-3'" = "TTN",
						 				 			 										 							 "Francisella FnCas12a (FnCpf1): 5'-YTN-3'" = "YTN"),
						 				 			 										 selected = "NGG"
						 				 			 	),
						 				 			 	actionLink("selectAll",  "Select All PAMs"),
						 				 			 	
						 				 			 	tags$br(),
						 				 			 	actionLink("selectNone", "De-select All PAMs"),
						 				 			 	
						 				 			 	p(paste0("If you select a PAM sequence which is a subset of another selected PAM sequence ", 
						 				 			 					 "(e.g., 'NGG' is a specific case of 'NRG'), MENTHU searches only for the more general case.")),
						 				 			 	
						 				 			 	tags$br(),
						 				 			 	textOutput("validpam"),
						 				 			 	
						 				 			 	#Input interface for custom pam sequences
						 				 			 	radioButtons("customCutOpt",
						 				 			 							 label    = "1a. Do you want to create a 'custom' PAM sequence?",
						 				 			 							 choices  = list("No"  = 0,
						 				 			 							 								 "Yes" = 1),
						 				 			 							 selected = 0, 
						 				 			 							 inline = TRUE),
						 				 			 	
						 				 			 	# When users use a custom input PAM sequence
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.customCutOpt == 1",
						 				 			 		textOutput("validmatchcustominputlength"),
						 				 			 		textOutput("validcustompam"),
						 				 			 		
						 				 			 		# tags$b(p("Input your nuclease PAM in the \"PAM_Sequence\" column. Ambiguous characters are allowed.")),
						 				 			 		# tags$b(p(paste0("Input the double-strand break position in relation to your PAM sequence in the \"DSB_Position\" column; ", 
						 				 			 		# 						"use negative values for DSBs upstream of the PAM, and positive values for downstream."))),
						 				 			 		# tags$b(p("Input the length of 5' overhang for sticky-cutting nucleases. Use '0' for blunt DSBs.")),
						 				 			 		# 
						 				 			 		# rHandsontableOutput("pamTable")
						 				 			 		
						 				 			 		# Area to input PAM sequences
						 				 			 		textAreaInput("customPamSeq",
						 				 			 									label = paste0("Input your nucleotide PAM sequence in the 5'-3' direction (e.g., 'NGG' for SpCas9).",
						 				 			 																 "You can enter multiple sequences by separating the PAMs with a comma or a space, e.g.",
						 				 			 																 "'NGG NRG NNVRYAC', etc.:"),
						 				 			 									value = "",
						 				 			 									placeholder = "Please input PAM sequence(s) in 5'-3' direction..."),
						 				 			 		textOutput("validcustomcutsites"),

						 				 			 		# Area to input cut patterns
						 				 			 		textAreaInput("cutSite",
						 				 			 									label = paste0("Please specify where your nuclease cuts in relation to your PAM ",
						 				 			 																 "(negative values for upstream). If you entered multiple PAM sequences,",
						 				 			 																 " please list the DSB locations in the ORDER YOU INPUT THE PAM SEQUENCES, e.g. '-3 -3 3', etc.:"),
						 				 			 									value = "",
						 				 			 									placeholder = "Input the DSB site(s) in relation to your PAM(s)..."),

						 				 			 		textAreaInput("overhang",
						 				 			 									label = paste0("Please specify the length of the overhangs resulting from the nuclease cut ",
						 				 			 																 "(use '0' (zero) for blunt cuts; a 5' overhange is assumed). If you entered multiple PAM sequences,",
						 				 			 																 " please list the overhang lengths in the ORDER YOU INPUT THE PAM SEQUENCES, e.g. '4 0 5', etc.:"),
						 				 			 									value = "",
						 				 			 									placeholder = "Input the overhang length(s)...")
						 				 			 	),
						 				 			 	
						 				 			 	tags$br(),
						 				 			 	
						 				 			 	#Decide whether to include TALEN targets
						 				 			 	radioButtons("talenOp",
						 				 			 							 label = "1b. Do you want to use TALEN targets?",
						 				 			 							 choices = list("No"  = 0,
						 				 			 							 							  "Yes" = 1),
						 				 			 							 selected = 0,
						 				 			 							 inline = TRUE
						 				 			 	),
						 				 			 	
						 				 			 	#For not using TALEN targets
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.talenOp == 0",
						 				 			 		p("You can use TALEN targets in addition to, or instead of, PAM targets.")
						 				 			 	),
						 				 			 	
						 				 			 	#Using TALEN targets
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.talenOp == 1",
						 				 			 		p(paste0("Specify TALEN design parameters. For an exact length, make \"Min\" and \"Max\" values identical.", 
						 				 			 						 "The \"Min\" value must be less than or equal to the \"Max\" value.")),
						 				 			 		textOutput("validtalen"),
						 				 			 		
						 				 			 		#Choose TALEN arm length by specifying min and max values
						 				 			 		tags$b("Choose the minimum and maximum TALEN arm length (arm length may range from 15-18 nucleotides): "),
						 				 			 		br(),
						 				 			 		div(style="display:inline-block",
						 				 			 				numericInput("armin",
						 				 			 										 label = "Min: ",
						 				 			 										 value = 15,
						 				 			 										 min   = 15,
						 				 			 										 max   = 18,
						 				 			 										 step  = 1,
						 				 			 										 width = "80px"
						 				 			 				)),
						 				 			 		
						 				 			 		div(style="display:inline-block",
						 				 			 				numericInput("armax",
						 				 			 										 label = "Max: ",
						 				 			 										 value = 15,
						 				 			 										 min   = 15,
						 				 			 										 max   = 18,
						 				 			 										 step  = 1,
						 				 			 										 width = "80px")
						 				 			 		),
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		tags$b("Choose the spacer length: "),
						 				 			 		tags$br(),
						 				 			 		div(style = "display:inline-block",
						 				 			 				radioButtons("spacer",
						 				 			 										 label = "",
						 				 			 										 choices = c("14 nts"       = 0,
						 				 			 										 						 "16 nts"       = 1,
						 				 			 										 						 "14 OR 16 nts" = 2),
						 				 			 										 inline = TRUE
						 				 			 				)
						 				 			 				
						 				 			 		)
						 				 			 	)
						 				 			 ),
						 				 			 
						 				 			 wellPanel(
						 				 			 	####Choose Input Type###############################################
						 				 			 	# Choose whether to use GenBank/RefSeq ID, Ensembl ID, or copy/paste input (used to correctly choose processing machinery)
						 				 			 	radioButtons("inputType",
						 				 			 							 label = "2. What type of sequence input do you want to use?",
						 				 			 							 choices = list("GenBank/RefSeq ID"        = 1,
						 				 			 							 							 "Ensembl ID"                = 3,
						 				 			 							 							 "Copy/paste gene sequence"  = 2
						 				 			 							 ),
						 				 			 							 selected = 1),
						 				 			 	
						 				 			 	####Input sequence#####################
						 				 			 	#If the user wants to use a genbank accession
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.inputType == 1",
						 				 			 		
						 				 			 		textAreaInput("genbankId",
						 				 			 									label       = paste0("Enter your GenBank or RefSeq ID here. ", 
						 				 			 																			 "For a full list of supported and not supported input types, ", 
						 				 			 																			 "please see the 'FAQs' section in the 'Instructions and FAQs' tab. "),
						 				 			 									value       = "",
						 				 			 									placeholder = "Enter GenBank gene ID here..."),
						 				 			 		
						 				 			 		# Validation check outputs
						 				 			 		textOutput("validgenbankid"),
						 				 			 		textOutput("genbankidexists")
						 				 			 		
						 				 			 	),
						 				 			 	
						 				 			 	# If user uses Ensembl accession
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.inputType == 3",
						 				 			 		textAreaInput("ensemblId",
						 				 			 									label       = paste0("Enter your Ensembl ID here. ",
						 				 			 																			 "We can analyze sequences from Ensembl transcript, exon, and protein IDs. ",
						 				 			 																			 "Please see the 'FAQs' section in the 'Instructions and FAQs' tab for more information."),
						 				 			 									value       = "",
						 				 			 									placeholder = "Enter Ensembl transcript, exon, or protein ID here..."),
						 				 			 		
						 				 			 		# Validation check outputs
						 				 			 		#textOutput("validensemblid"),
						 				 			 		textOutput("ensemblidexists")
						 				 			 		
						 				 			 	),
						 				 			 	
						 				 			 	#If user wants to copy/paste sequence
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.inputType == 2",
						 				 			 		
						 				 			 		#Text area to copy/paste gene sequence
						 				 			 		textAreaInput("geneSeq",
						 				 			 									label       = "Please input your genomic sequence of interest (sequences should be between 80 and 5000 nucleotides in length.)",
						 				 			 									value       = "",
						 				 			 									placeholder = "Paste gene sequence here..."),
						 				 			 		
						 				 			 		#Check to make sure gene sequence is valid
						 				 			 		textOutput("validgeneseq"), 
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		#Choose to specify exon inputs by hand
						 				 			 		radioButtons("pasteExonType", 
						 				 			 								 label    = "Does your pasted sequence have multiple exons, and do you wish to find a target within those exons?",
						 				 			 								 choices  = list("No"  = 0, 
						 				 			 								 								"Yes" = 1),
						 				 			 								 selected = 0, 
						 				 			 								 inline   = TRUE),
						 				 			 		
						 				 			 		#If the user wants to enter exon information
						 				 			 		conditionalPanel(
						 				 			 			condition = "input.pasteExonType == 1",
						 				 			 			p(paste0("Input exon location information here. Only include information about exons you wish to target.", 
						 				 			 							 "If your first exon of interest starts on nucleotide 63 and ends on nucleotide 145 of your sequence,",
						 				 			 							 " type '63' into 'exonStart' and '145' into 'exonEnd'.")),
						 				 			 			p(paste0("If you need to add more exons, right-click or command-click and select 'Insert Row'.", 
						 				 			 							 " You do not need to remove extra rows.")),
						 				 			 			rHandsontableOutput("exonInfo"),
						 				 			 			textOutput("validexoninfo")
						 				 			 		)
						 				 			 	)
						 				 			 ),
						 				 			 
						 				 			 ############Exon Options#########################
						 				 			 conditionalPanel(
						 				 			 	condition = "input.inputType == 1 | input.inputType == 3",
						 				 			 	
						 				 			 	wellPanel(
						 				 			 		#Choose whether to include first exon in future determinations
						 				 			 		radioButtons("firstExon",
						 				 			 								 label    = "2a. Do you want to find targets in the first exon? (Not recommended for gene knockouts.)",
						 				 			 								 choices  = list("No"  = 0,
						 				 			 								 								"Yes" = 1),
						 				 			 								 selected = 0,
						 				 			 								 inline   = TRUE),
						 				 			 		
						 				 			 		#Choose option for specifying exon input
						 				 			 		p(tags$b("Which exon(s) do you want to target?")),
						 				 			 		radioButtons("exonTargetType",
						 				 			 								 label   = "Do you want to: ",
						 				 			 								 choices = list("Search all exons - Please note if 2a = 'No', the first exon will not be searched"           = 0,
						 				 			 								 							 "Search for a target within a specified percentage of exons of the beginning of the sequence" = 1,
						 				 			 								 							 "Search for a target within a specified percentage of exons of the end of the sequence"       = 2,
						 				 			 								 							 "Provide a list of exons to target"                                                           = 3),
						 				 			 								 selected = 0),
						 				 			 		
						 				 			 		#When working from the beginning of the sequence
						 				 			 		conditionalPanel(
						 				 			 			condition = "input.exonTargetType == 1",
						 				 			 			p(paste0("Specify a percentage of exons, starting from the beginning of the sequence, ", 
						 				 			 							 "to search (this includes the first exon if you chose 'yes' in 3, and excludes it if you chose 'no'.) ", 
						 				 			 							 "For gene knockouts, it is recommended that you use 30%. ",
						 				 			 							 "For example, if your gene has 10 exons and you choose not to look for targets in the first exon, ", 
						 				 			 							 "a value of '30%' will search exons 2, 3, and 4 for target sites. ",
						 				 			 							 "If you chose to look for targets in the first exon, a value of 30% will search exons 1, 2, and 3. ")),
						 				 			 			p("You can look for targets in every exon (100%) by changing the value to '100'. Minimum value is 1%."),
						 				 			 			
						 				 			 			numericInput("exonBegPercentage",
						 				 			 									 label = "",
						 				 			 									 min   = 1,
						 				 			 									 max   = 100,
						 				 			 									 value = 30)
						 				 			 		),
						 				 			 		
						 				 			 		# When working from the end of the sequence
						 				 			 		conditionalPanel(
						 				 			 			condition = "input.exonTargetType == 2",
						 				 			 			p("Specify a percentage of exons, starting from the end of the sequence, to search. "),
						 				 			 			
						 				 			 			numericInput("exonEndPercentage",
						 				 			 									 label = "",
						 				 			 									 value = 30,
						 				 			 									 min   = 1,
						 				 			 									 max   = 100)
						 				 			 		),
						 				 			 		
						 				 			 		#When listing exons to target
						 				 			 		conditionalPanel(
						 				 			 			condition = "input.exonTargetType == 3",
						 				 			 			textAreaInput("exonTargetList",
						 				 			 										label = paste0("Specify the exon(s) you wish to target. ",
						 				 			 																	 "Multiple exons should be separated with a comma ",
						 				 			 																	 "(e.g., to select exon 3 and 5, type '3,5'.) ",
						 				 			 																	 "A range of exons can be specified by a dash ",
						 				 			 																	 "(e.g., to target all exons between exon 5 and exon 10, type '5-10'). ",
						 				 			 																	 "You can specify multiple ranges and/or multiple exons using commas ",
						 				 			 																	 "(e.g., to select exon 3, exon 5 through 7, exon 9, ",
						 				 			 																	 "and exon 11 through 13, type '3,5-7,9,11-13'.)"), 
						 				 			 										placeholder = "Enter exons here...",
						 				 			 										value       = ""),
						 				 			 			textOutput("validexonlist")
						 				 			 		)
						 				 			 	)#,
						 				 			 	
						 				 			 	#wellPanel(
						 				 			 	#	
						 				 			 	# Allow user to choose if context consideration can run over into exon region
						 				 			 	#	radioButtons("contextWiggleType",
						 				 			 	#							 label    = paste0("2c. Do you want to include target sites where the double-strand break site will occur within an exon, ",
						 				 			 	#							 									"but the contextual information used to calculate the MENTHU score may include intronic sequences?"),
						 				 			 	#							 choices  = c("No" = 0,
						 				 			 	#							 						 "Yes" = 1),
						 				 			 	#							 selected = 0,
						 				 			 	#							 inline   = TRUE
						 				 			 	#	),
						 				 			 	
						 				 			 	# Allow user to choose if gRNA can run over into exon region
						 				 			 	#	radioButtons("gRNAWiggleType",
						 				 			 	#							 label    = paste0("2b. Do you want to include target sites where the double-strand break ",
						 				 			 	#							 							     "site will occur within an exon, but the gRNA may run over the exon boundary into an intron?"),
						 				 			 	#							 choices  = c("No"  = 0,
						 				 			 	#							 						  "Yes" = 1),
						 				 			 	#							 selected = 0,
						 				 			 	#							 inline   = TRUE
						 				 			 	#	)
						 				 			 	#	
						 				 			 	#	
						 				 			 	#	
						 				 			 	#)
						 				 			 ),
						 				 			 
						 				 			 wellPanel(	
						 				 			 	
						 				 			 	#Submit panel for GenBank submissions
						 				 			 	conditionalPanel(
						 				 			 		# If GenBank/RefSeq
						 				 			 		condition = "input.inputType == 1", 
						 				 			 		
						 				 			 		#Submit button
						 				 			 		actionButton("genbankSubmit", "Submit"), 
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		p("Your results may take a few minutes to calculate. Please do not close this web page until your calculation is finished."),
						 				 			 		
						 				 			 		# Generate the UI for the download button
						 				 			 		uiOutput("downOutGB"),
						 				 			 		tags$br(), 
						 				 			 		
						 				 			 		# Generate button to bookmark inputs
						 				 			 		bookmarkButton(label = "Bookmark Session Inputs",
						 				 			 									 title = "Click here to generate a URL that can be copy/pasted into a web browser to easily reproduce your analysis",
						 				 			 									 id    = "bookmarkGB"),
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		#Indicates if there is a problem with the supplied ID
						 				 			 		textOutput('genbankIdOutcome'), 
						 				 			 		tags$style("#genbankIdOutcome{color: red;}"),
						 				 			 		
						 				 			 		# Output information regarding the number of hits and misses
						 				 			 		#uiOutput('genbankhits'),
						 				 			 		#tags$br(),
						 				 			 		
						 				 			 		#Output results
						 				 			 		conditionalPanel(
						 				 			 			condition = "output.filtOpsGB",
						 				 			 			tags$br(),
						 				 			 			p("Filter Options: "),
						 				 			 			checkboxInput("t7OptGB",     "T7-compatible gRNAs",                       value = FALSE),
						 				 			 			checkboxInput("thresholdGB", "Recommended sites (>=1.5 score threshold)", value = FALSE)
						 				 			 		),
						 				 			 		
						 				 			 		uiOutput('genbankResults')
						 				 			 	),
						 				 			 	
						 				 			 	#Submit panel for Ensembl submissions
						 				 			 	conditionalPanel(
						 				 			 		# If Ensembl
						 				 			 		condition = "input.inputType == 3",
						 				 			 		
						 				 			 		# Submit button
						 				 			 		actionButton("ensemblSubmit", "Submit"), 
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		p("Your results may take a few minutes to calculate. Please do not close this web page until your calculation is finished."),
						 				 			 		
						 				 			 		# Generate UI for download button
						 				 			 		uiOutput("downOutEns"),
						 				 			 		tags$br(), 
						 				 			 		
						 				 			 		# Generate button to bookmark inputs
						 				 			 		bookmarkButton(label = "Bookmark Session Inputs",
						 				 			 									 title = "Click here to generate a URL that can be copy/pasted into a web browser to easily reproduce your analysis",
						 				 			 									 id    = "bookmarkE"),
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		# Indicates if there is a problem with the supplied ID or connection to Ensembl
						 				 			 		#textOutput('ensemblIdOutcome'), 
						 				 			 		#tags$style("#ensemblIdOutcome{color: red;}"),
						 				 			 		textOutput('ensemblUp'),
						 				 			 		tags$style("#ensemblUp{color: red;"),
						 				 			 		
						 				 			 		#uiOutput('ensemblhits'),
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		#Output results
						 				 			 		conditionalPanel(
						 				 			 			condition = "output.filtOpsE",
						 				 			 			tags$br(),
						 				 			 			p(tags$b("Filter Options: ")),
						 				 			 			checkboxInput("t7OptE",     "T7-compatible gRNAs",                       value = FALSE),
						 				 			 			checkboxInput("thresholdE", "Recommended sites (>=1.5 score threshold)", value = FALSE)
						 				 			 		),
						 				 			 		
						 				 			 		uiOutput('ensemblResults')
						 				 			 	),
						 				 			 	
						 				 			 	conditionalPanel(
						 				 			 		condition = "input.inputType == 2",
						 				 			 		actionButton("geneSeqSubmit", "Submit"),
						 				 			 		
						 				 			 		tags$br(),
						 				 			 		p("Your results may take a few minutes to appear. Please do not close this web page until your calculation is finished."),
						 				 			 		uiOutput("downOutGS"),
						 				 			 		tags$br(), 
						 				 			 		bookmarkButton(label = "Bookmark Session Inputs",
						 				 			 									 title = "Click here to generate a URL that can be copy/pasted into a web browser to easily reproduce your analysis",
						 				 			 									 id    = "bookmarkGS"),
						 				 			 		tags$br(),
						 				 			 		
						 				 			 		textOutput('geneseqerrors'),
						 				 			 		
						 				 			 		#Output results
						 				 			 		#uiOutput('geneseqhits'),
						 				 			 		#br(),
						 				 			 		conditionalPanel(
						 				 			 			condition = "output.filtOpsGS",
						 				 			 			tags$br(),
						 				 			 			p("Filter Options: "),
						 				 			 			checkboxInput("t7OptGS",     "T7-compatible gRNAs",                       value = FALSE),
						 				 			 			checkboxInput("thresholdGS", "Recommended sites (>=1.5 score threshold)", value = FALSE)
						 				 			 		),
						 				 			 		
						 				 			 		uiOutput('geneSeqResults')
						 				 			 		#tags$head(tags$style("#geneSeqResults table {background-color: white; }", media = "screen", type = "text/css"))
						 				 			 		
						 				 			 	)
						 				 			 )
						 				 )
						 ),
						 
						 ##########Pre-COMPUTED GENES TAB#################################
						 # This is being moved to its own app
						 
						 ##########TOOLS AND DOWNLOADS TAB#################################
						 
						 tabPanel(
						 	tags$div("Tools and Downloads"),
						 	titlePanel(""),
						 	
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 												 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 												 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 												 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://github.com/Dobbs-Lab/MENTHU", 
						 												 target = "_blank", tags$img(src = "GitHub_Logo.png",                   width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://hub.docker.com/r/cmmann/menthu", 
						 												 target = "_blank", tags$img(src = "Docker_Logo.png",                   width = "100%"))))
						 		
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		h3("Download MENTHU"),
						 		tags$p(HTML(paste0("A standalone version of this code can be downloaded from our ", 
						 											 tags$a(href = "https://github.com/Dobbs-Lab/MENTHU", target = "_blank", "GitHub repository"),
						 											 "."))),
						 		tags$p(HTML(paste0("There are extensive installation/usage instructions available in the GitHub ", 
						 											 tags$a(href = "https://github.com/Dobbs-Lab/MENTHU#how-to-run-menthu-locally", target = "_blank", "README"), 
						 											 " file."))),
						 		tags$p("You can clone the repository with the following git command:"),
						 		tags$p(tags$code("git clone https://github.com/Dobbs-Lab/MENTHU.git"), style = "text-align:center;"),
						 		tags$p(HTML(paste0("MENTHU is available under the GNU General Public License v3 (GPL 3.0). You can read the license ",
						 											 tags$a(href = "https://github.com/Dobbs-Lab/MENTHU/blob/master/LICENSE", target = "_blank", "here"),
						 											 "."))),
						 		tags$p(HTML(paste0("The MENTHU R code is provided as-is; please be aware that you modify the code at your own risk. ",
						 											 "We are unable to provide technical support for modified versions.")))
						 	),
						 		
						 		wellPanel(
						 			h3("Run MENTHU Locally"),
						 		tags$p(HTML(paste0("If you have R installed on your system, you can also follow the instructions ",
						 											 tags$a(href = "https://github.com/Dobbs-Lab/MENTHU#3-run-menthu-locally", target = "_blank", "here"),
						 											 " to easily run the MENTHU RShiny app from R, without dealing with Git."))),
						 		p("MENTHU is also available as a Docker container image. You can clone the Docker image using the following command:"),
						 		tags$p(tags$code("sudo docker pull cmmann/menthu"), style = "text-align:center;")
						 		
						 	))
						 	
						 ),
						 
						 ##########FUNDING Tab#############################################
						 tabPanel(
						 	tags$div("Funding"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		#Sidebar panel with links
						 		tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 												 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 												 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 												 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 												 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 												 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 												 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://github.com/Dobbs-Lab/MENTHU", 
						 												 target = "_blank", tags$img(src = "GitHub_Logo.png",                   width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://hub.docker.com/r/cmmann/menthu", 
						 												 target = "_blank", tags$img(src = "Docker_Logo.png",                   width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://dill-picl.org", 
						 												 target = "_blank", tags$img(src = "lawlab_web_wiki_header.png",        width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.nih.gov/", 
						 												 target = "_blank", tags$img(src = "nihlogo.png",                       width = "100%"))))
						 		
						 	)),
						 	
						 	column(9, wellPanel(
						 		includeHTML("www/funding.html")
						 	))
						 ),
						 
						 
						 ##########HOW TO CITE Tab#########################################
						 tabPanel(
						 	tags$div("How to Cite"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 												 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 												 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 												 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://github.com/Dobbs-Lab/MENTHU", 
						 												 target = "_blank", tags$img(src = "GitHub_Logo.png",                   width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://hub.docker.com/r/cmmann/menthu", 
						 												 target = "_blank", tags$img(src = "Docker_Logo.png",                   width = "100%"))))
						 		
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		includeHTML("www/citation.html")
						 	))
						 	
						 ),
						 
						 ##########CONTACT US Tab##########################################
						 tabPanel(
						 	tags$div("Report Bugs or Contact Us"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 												 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 												 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 												 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://github.com/Dobbs-Lab/MENTHU", 
						 												 target = "_blank", tags$img(src = "GitHub_Logo.png",                   width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://hub.docker.com/r/cmmann/menthu", 
						 												 target = "_blank", tags$img(src = "Docker_Logo.png",                   width = "100%"))))
						 		
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		p("Please use the form below, or email us directly at help@genesculpt.org, to report issues and request support.")
						 	),
						 	
						 	tags$iframe(id           = "googleform", 
						 							src          = paste0("https://docs.google.com/forms/d/e/1FAIpQLSeq9aDRj6EOCskBwPsA2PFQ2LsKxT4v85-", 
						 																		"rGTlYQOk0n8X2Gw/viewform?usp=pp_url&entry.358268393&entry.1646278736=MENTHU&entry.", 
						 																		"1934309806&entry.565411344&entry.754537383&entry.826100992"),
						 							width        = 760,
						 							height       = 2000,
						 							frameborder  = 0,
						 							marginheight = 0)
						 	)
						 ),
						 
						 ######STATUS and CHANGELOG Tab####################################
						 tabPanel(
						 	tags$div("Change Log"),
						 	titlePanel(""),
						 	
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href   = "http://genesculpt.org/gss/", 
						 												 target = "_blank", tags$img(src = "GSS logo small.png",                width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.iastate.edu/",   
						 												 target = "_blank", tags$img(src = "isu-logo-alt.png",                  width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.mayoclinic.org", 
						 												 target = "_blank", tags$img(src = "MC_stack_4c_DAC.png",               width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://github.com/Dobbs-Lab/MENTHU", 
						 												 target = "_blank", tags$img(src = "GitHub_Logo.png",                   width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href   = "https://hub.docker.com/r/cmmann/menthu", 
						 												 target = "_blank", tags$img(src = "Docker_Logo.png",                   width = "100%"))))
						 		
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		includeHTML("www/changelog.html")
						 	))
						 )
	)}
)

#