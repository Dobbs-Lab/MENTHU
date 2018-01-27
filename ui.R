#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(rhandsontable)
library(shinyIncubator)

shinyUI(
	
	####Creates the navbar set-up####
	navbarPage(id = 'mainPage',
						 windowTitle = "MENTHU",
						 
						 #Stylesheet
						 theme = "ogtheme.css", 
						 
						 #Page title box
						 tags$div("MENTHU v1.1.0", 
						 				 style = "color:white"),
						 
						 ########ABOUT TAB#################################################
						 tabPanel(tags$div("About", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 				 											 target = "_blank", 
						 				 											 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 				 )),
						 				 
						 				 #Text area in center of page
						 				 column(9, wellPanel(
						 				 	includeHTML("www/menthuAbout.html")
						 				 ))
						 				 
						 ),
						 
						 ##########INSTRUCTIONS############################################
						 tabPanel(tags$div("Instructions and FAQs", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 				 	tags$br(),
						 				 	tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 				 											 target = "_blank", 
						 				 											 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 				 )),
						 				 
						 				 #Text area in center of page
						 				 column(9, wellPanel(
						 				 	
						 				 	#Display page instructions
						 				 	includeHTML("www/menthuInstructions.html")
						 				 ))
						 				 
						 ),
						 
						 ##########Calculate MENTHU TAB###################################
						 tabPanel(id = "single",
						 				 tags$div("Calculate MENTHU Score", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 ##Sidebar############################################################
						 				 #Adds a sidebar for users to pre-populate fields with an example, and reset the form
						 				 column(2, wellPanel(
						 				 	
						 				 	#GenBank Example
						 				 	actionLink("exampleGenBank",
						 				 						 label = "GenBank Gene ID Example"),
						 				 	p(""),
						 				 	
						 				 	#Cut/Paste cDNA example; input$example
						 				 	actionLink("exampleGeneSeq", 
						 				 						 label = "Pasted Sequence Example"),
						 				 	
						 				 	p(""),
						 				 	
						 				 	#Reset Button; input$reset
						 				 	actionLink("reset", 
						 				 						 label = "Reset Form")
						 				 )),
						 				 
						 				 
						 				 ####Main Bar#########################################################
						 				 #Main panel for entering information and submitting job
						 				 column(9, wellPanel(
						 				 	
						 				 	####Choose PAM sequence#############################################
						 				 	checkboxGroupInput("casType",
						 				 										 label = "1. Select the PAM sequence(s) you wish to target:",
						 				 										 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = "NGG",
						 				 										 							 "S. pyogenes SpCas9: 5'-NRG-3'" = "NRG",
						 				 										 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = "NNNRRT",
						 				 										 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = "NNGRRT",
						 				 										 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = "NNGTGA",
						 				 										 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = "NNAGAAW",
						 				 										 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = "NNNVRYAC",
						 				 										 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = "NNNNGMTT"),
						 				 										 #CPF1 is not currently supported, but may be in future releases
						 				 										 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = "TTTN",
						 				 										 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = "TTTV",
						 				 										 #"Francisella FnCpf1: 5'-TTN-3'" = "TTN",
						 				 										 #"Francisella FnCpf1: 5'-YTN-3'" = "YTN"),
						 				 										 selected = "NGG"
						 				 	),

						 				 	
						 				 	tags$br(),
						 				 	actionLink("selectAll", "Select All PAMs"),
						 				 	tags$br(),
						 				 	actionLink("selectNone", "De-select All PAMs"),
						 				 	tags$br(),
						 				 	textOutput("validpam"),
						 				 	
						 				 	
						 				 	radioButtons("customCutOpt",
						 				 							 label    = "1a. Do you want to create a 'custom' PAM sequence?",
						 				 							 choices  = list("No" = 0,
						 				 							 								"Yes" = 1),
						 				 							 selected = 0, 
						 				 							 inline = TRUE),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.customCutOpt == 1",
												textOutput("validmatchcustominputlength"),
												textOutput("validcustompam"),
						 				 		textAreaInput("customPamSeq",
						 				 									label = "Input your nucleotide PAM sequence in the 5'-3' direction (e.g., 'NGG' for SpCas9). You can enter multiple sequences by separating the PAMs with a comma or a space, e.g. 'NGG NRG NNVRYAC', etc.:", 
						 				 									#\nIf your sequence allows for ambiguity, please use the IUPAC one-letter codes: R = A or G; Y = C or T; S = G or C; W = A or T; K = G or T; M = A or C; B = not A; D = not C; H = not G; V = not T; N = any base.",
						 				 									value = "",
						 				 									placeholder = "Please input PAM sequence(s) in 5'-3' direction..."),
												textOutput("validcustomcutsites"),
						 				 		textAreaInput("cutSite",
						 				 									label = "Please specify where your nuclease cuts in relation to your PAM (negative values for upstream). If you entered multiple PAM sequences, please list the DSB locations in the ORDER YOU INPUT THE PAM SEQUENCES, e.g. '-3 -3 3', etc.:",
#If the double-strand break site (DSB) is 5' to the PAM, count the number of bases between the DSB site and the FIRST base of your PAM. For example, SpCas9 generally induces a DSB three bases upstream of 'NGG', so the input value would be -3. If the cut is 3' to the PAM, count the number of bases between the LAST base of your PAM and the DSB site.",
						 				 									value = "",
						 				 									placeholder = "Input the DSB site(s) in relation to your PAM(s)...")
						 				 	),
						 				 	
						 				 	#Decide whether to include TALEN targets
						 				 	radioButtons("talenOp",
						 				 							 label = "1b. Do you want to use TALEN targets?",
						 				 							 choices = list(
						 				 							 	"No" = 0,
						 				 							 	"Yes" = 1
						 				 							 ),
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
						 				 						 "The \"Min\" value must be less than or equal to the \"Max\" value. Please note that to shorten computation time, we currently only search for TALENs in the forward/sense orientation.")),
						 				 		textOutput("validtalen"),
						 				 		
						 				 		#Choose TALEN arm length by specifying min and max values
						 				 		tags$b("Choose the minimum and maximum TALEN arm length (arm length may range from 15-18 nucleotides): "),
						 				 		tags$br(),
						 				 		div(style="display:inline-block",
						 				 				numericInput("armin",
						 				 										 label = "Min: ",
						 				 										 value = 15,
						 				 										 min = 15,
						 				 										 max = 18,
						 				 										 step = 1,
						 				 										 width = "80px"
						 				 				)),
						 				 		
						 				 		div(style="display:inline-block",
						 				 				numericInput("armax",
						 				 										 label = "Max: ",
						 				 										 value = 15,
						 				 										 min = 15,
						 				 										 max = 18,
						 				 										 step = 1,
						 				 										 width = "80px")
						 				 		),
						 				 		tags$br(),
						 				 		
						 				 		#Choose spacer length by specifying min and max values
						 				 		tags$b("Choose the minimum and maximum spacer length (spacer length may range from 14-16 nucleotides): "),
						 				 		tags$br(),
						 				 		div(style="display:inline-block",
						 				 				numericInput("spamin",
						 				 										 label = "Min: ",
						 				 										 value = 15,
						 				 										 min = 14,
						 				 										 max = 16,
						 				 										 step = 1,
						 				 										 width = "80px"
						 				 				)),
						 				 		
						 				 		div(style="display:inline-block",
						 				 				numericInput("spamax",
						 				 										 label = "Max: ",
						 				 										 value = 15,
						 				 										 min = 14,
						 				 										 max = 16,
						 				 										 step = 1,
						 				 										 width = "80px")
						 				 		)
						 				 	)
						 				 ),
						 				 
						 				 wellPanel(
						 				 	####Choose Input Type###############################################
						 				 	radioButtons("inputType",
						 				 							 label = "2. What type of sequence input do you want to use?",
						 				 							 choices = list("GenBank Gene ID"  = 1,
						 				 							 							 "Copy/paste gene sequence" = 2
						 				 							 ),
						 				 							 selected = 1),
						 				 	
						 				 	####Input sequence#####################
						 				 	#If the user wants to use a genbank accession
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 1",
						 				 		
						 				 		textAreaInput("genbankId",
						 				 									label = "",
						 				 									value = NULL,
						 				 									placeholder = "Paste GenBank gene ID here..."),
						 				 		textOutput("validgenbankid"),
						 				 		textOutput("genbankidexists")
						 				 		
						 				 	),
						 				 	
						 				 	#If user wants to copy/paste sequence
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 2",
						 				 		
						 				 		#Text area to copy/paste gene sequence
						 				 		textAreaInput("geneSeq",
						 				 									label = "",
						 				 									value = "",
						 				 									placeholder = "Paste gene sequence here..."),
						 				 		
						 				 		#Check to make sure gene sequence is valid
						 				 		textOutput("validgeneseq"), 
						 				 		tags$br(),
						 				 		
						 				 		#Choose to specify exon inputs by hand
						 				 		radioButtons("pasteExonType", 
						 				 								 label = "Does your pasted sequence have multiple exons, and do you wish to find a target within those exons?",
						 				 								 choices = list("No" = 0, "Yes" = 1),
						 				 								 selected = 0, 
						 				 								 inline = TRUE),
						 				 		
						 				 		#If the user wants to enter exon information
						 				 		conditionalPanel(
						 				 			condition = "input.pasteExonType == 1",
						 				 			p(paste0("Input exon location information here. Only include information about exons you wish to target.", 
						 				 							 "If your first exon of interest starts on nucleotide 63 and ends on nucleotide 145 of your sequence,",
						 				 							 " type '63' into 'exonStart' and '145' into 'exonEnd'.")),
						 				 			p(paste0("If you need to add more exons, right-click or command-click and select 'Insert Row'.", 
						 				 							 " You do not need to remove extra rows.")),
						 				 			#uiOutput("inputTable"),
						 				 			rHandsontableOutput("exonInfo"),
						 				 			textOutput("validexoninfo")
						 				 		)
						 				 	)
						 				 ),
						 				 
						 				 
						 				 ############Exon Options#########################
						 				 conditionalPanel(
						 				 	condition = "input.inputType == 1",
						 				 	
						 				 	wellPanel(
						 				 		#Choose whether to include first exon in future determinations
						 				 		radioButtons("firstExon",
						 				 								 label = "2a. Do you want to find targets in the first exon? (Not recommended for gene knockouts.)",
						 				 								 choices = list("No" = 0,
						 				 								 							 "Yes" = 1),
						 				 								 selected = 0,
						 				 								 inline = TRUE),
						 				 		
						 				 		#Choose option for specifying exon input
						 				 		p(tags$b("Which exon(s) do you want to target?")),
						 				 		radioButtons("exonTargetType",
						 				 								 label = "Do you want to: ",
						 				 								 choices = list("Search all exons" = 0,
						 				 								 							 "Search for a target within a specified percentage of exons of the beginning of the sequence" = 1,
						 				 								 							 "Search for a target within a specified percentage of exons of the end of the sequence" = 2,
						 				 								 							 "Provide a list of exons to target" = 3),
						 				 								 selected = 0),
						 				 		
						 				 		#When working from the beginning of the sequence
						 				 		conditionalPanel(
						 				 			condition = "input.exonTargetType == 1",
						 				 			p(paste0("Specify a percentage of exons, starting from the beginning of the sequence,", 
						 				 							 " to search (this includes the first exon if you chose 'yes' in 3, and excludes it if you chose 'no'.)", 
						 				 							 " For gene knockouts, it is recommended that you use 30%. ",
						 				 							 "For example, if your gene has 10 exons and you choose not to look for targets in the first exon,", 
						 				 							 " a value of '30%' will search exons 2, 3, and 4 for target sites. ",
						 				 							 "If you chose to look for targets in the first exon, a value of 30% will search exons 1, 2, and 3.")),
						 				 			p("You can look for targets in every exon (100%) by changing the value to '100'. Minimum value is 1%."),
						 				 			
						 				 			numericInput("exonBegPercentage",
						 				 									 label = "",
						 				 									 min = 1,
						 				 									 max = 100,
						 				 									 value = 30)
						 				 		),
						 				 		
						 				 		#When working from the end of the sequence
						 				 		conditionalPanel(
						 				 			condition = "input.exonTargetType == 2",
						 				 			p("Specify a percentage of exons, starting from the end of the sequence, to search. "),
						 				 			
						 				 			numericInput("exonEndPercentage",
						 				 									 label = "",
						 				 									 value = 30,
						 				 									 min = 1,
						 				 									 max = 100)
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
						 				 										placeholder = "Enter exons here...")
						 				 		)
						 				 	)
						 				 ),
						 				 
						 				 ###########THRESHOLD#############################
						 				 wellPanel(
						 				 	#Choose value for MENTHU display threshold
						 				 	numericInput("threshold",
						 				 							 label = paste0("3. Choose minimum score for reporting (>50 is a strong score, ",
						 				 							 							 ">40 is moderate; we do not recommend using sites <40.): "),
						 				 							 value = 40,
						 				 							 min = 0,
						 				 							 max = 100
						 				 	),
						 				 	
						 				 	#Validate threshold
						 				 	textOutput("validthreshold")
						 				 ),
						 				 
						 				 wellPanel(	
						 				 	
						 				 	#Submit panel for GenBank subissions
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 1",
						 				 		actionButton("genbankSubmit", "Submit"), #Submit button
						 				 		tags$br(),
						 				 		
						 				 		uiOutput("downOutGenbank"), #Download button; should only display when there are results to download
						 				 		tags$br(),
						 				 		
						 				 		textOutput('genbankIdOutcome'), #Indicates if there is a problem with the supplied s
						 				 		tags$style("#genbankIdOutcome{color: red;}"),
						 				 		DT::dataTableOutput('genbankResults')
						 				 		
						 				 	),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 2",
						 				 		actionButton("geneSeqSubmit", "Submit"),
						 				 		tags$br(),
						 				 		
						 				 		uiOutput("downOut"),
						 				 		tags$br(),
						 				 		
						 				 		DT::dataTableOutput('geneSeqResults')
						 				 		#tags$head(tags$style("#geneSeqResults table {background-color: white; }", media = "screen", type = "text/css"))
						 				 	),
						 				 	p("Your results may take a few minutes to appear. Please do not close this web page until your calculation is finished.")
						 				 	
						 				 )
						 				 )
						 ),
						 
						 
						 ##########HOW TO CITE Tab#########################################
						 tabPanel(
						 	tags$div("How to Cite", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", 
						 												 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		p("Manuscript is in prep; citation will be available shortly.")
						 	))
						 	
						 ),
						 
						 ##########CONTACT US Tab##########################################
						 tabPanel(
						 	tags$div("Report Bugs or Request Help", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", 
						 												 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		p("Please use the form below, or email us directly at GeneSculptSuite@gmail.com, to report issues and request support.")
						 	),
						 	
						 	tags$iframe(id = "googleform", 
						 							src = "https://docs.google.com/forms/d/e/1FAIpQLSeq9aDRj6EOCskBwPsA2PFQ2LsKxT4v85-rGTlYQOk0n8X2Gw/viewform?usp=pp_url&entry.358268393&entry.1646278736=MENTHU&entry.1934309806&entry.565411344&entry.754537383&entry.826100992",
						 							width = 760,
						 							height = 2000,
						 							frameborder = 0,
						 							marginheight = 0)
						 	)
						 	
						 ),
						 ##########TOOLS AND DOWNLOADS TAB#################################
						 
						 tabPanel(
						 	tags$div("Tools and Downloads", style = "color:white"),
						 	titlePanel(""),
						 	
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", 
						 												 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		p("A standalone version of this code may be downloaded from", tags$a(href = "https://github.com/Dobbs-Lab/MENTHU", target = "_blank", " GitHub."), " The R code is provided as-is, and may not be used in commercial applications. Please be aware that you modify the code at your own risk; we are unable to provide support for modified versions.")
						 	))
						 	
						 ),
						 
						 ##########FUNDING Tab#############################################
						 tabPanel(
						 	tags$div("Funding", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		#Sidebar panel with links
						 		tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", 
						 												 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.nih.gov/", target = "_blank", tags$img(src = "nihlogo.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://dill-picl.org", target = "_blank", tags$img(src = "lawlab_web_wiki_header.png", width = "100%"))))
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		tags$p("This webtool was created and is maintained by funding through NIH R24 OD020166 and is a joint effort by Iowa State University and The Mayo Clinic."),
						 		tags$p("This server is generously hosted by the", a(href = "https://dill-picl.org/", target = "_blank", 'Lawrence-Dill Plant Informatics and Computation (Dill-PICL) Lab'), "at Iowa State University.")
						 	))
						 ),
						 
						 
						 ######STATUS and CHANGELOG Tab####################################
						 tabPanel(
						 	tags$div("Change Log", style = "color:white"),
						 	titlePanel(""),
						 	
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		tags$div(tags$span(a(href = "http://ll-g2f.gdcb.iastate.edu/gss/", target = "_blank", tags$img(src = "GSS logo small.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isu-logo-alt.png",     width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "MC_stack_4c_DAC.png", width = "100%")))),
						 		tags$br(),
						 		tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", 
						 												 target = "_blank", 
						 												 tags$img(src = "genome-writers-guild-logo_DAC.png", width = "100%"))))
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		includeHTML("www/changelog.html")
						 	))
						 )
	))
#