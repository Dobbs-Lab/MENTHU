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
						 tags$div("MENTHU v1.0.1", 
						 				 style = "color:white"),
						 
						 ########ABOUT TAB#################################################
						 tabPanel(tags$div("About", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href = "https://www.iastate.edu/",   
						 				 											 target = "_blank", 
						 				 											 tags$img(src = "isulogo.jpg",     
						 				 											 				 width = "90%")))),
						 				 	p(""),
						 				 	tags$div(tags$span(a(href = "https://www.mayoclinic.org", 
						 				 											 target = "_blank", 
						 				 											 tags$img(src = "mayoclinic.jpeg", 
						 				 											 				 width = "50%"))))
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
						 				 	tags$div(tags$span(a(href = "https://www.iastate.edu/",   
						 				 											 target = "_blank", 
						 				 											 tags$img(src = "isulogo.jpg",    
						 				 											 				 width = "90%")))),
						 				 	p(""),
						 				 	tags$div(tags$span(a(href = "https://www.mayoclinic.org", 
						 				 											 target = "_blank", 
						 				 											 tags$img(src = "mayoclinic.jpeg", 
						 				 											 				 width = "50%"))))
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
						 				 						 label = "Copy/paste Gene Sequence Example"),
						 				 	
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
						 				 	
						 				 	p(""),
						 				 	actionLink("selectAll", "Select All PAMs"),
						 				 	p(""),
						 				 	actionLink("selectNone", "De-select All PAMs"),
						 				 	p(""),
						 				 	textOutput("validpam"),
						 				 	
						 				 	#Decide whether to include TALEN targets
						 				 	radioButtons("talenOp",
						 				 							 label = "1a. Do you want to use TALEN targets?",
						 				 							 choices = list(
						 				 							 	"No" = 0,
						 				 							 	"Yes" = 1
						 				 							 ),
						 				 							 selected = 0,
						 				 							 inline = TRUE
						 				 	),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.talenOp == 0",
						 				 		p("You can use TALEN targets in addition to, or instead of, PAM targets.")
						 				 	),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.talenOp == 1",
						 				 		p("Specify TALEN design parameters. For an exact length, make \"Min\" and \"Max\" values identical. The \"Min\" value must be less than or equal to the \"Max\" value."),
						 				 		#Choose TALEN arm length by specifying min and max values
						 				 		tags$b("Choose the minimum and maximum TALEN arm length: "),
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
						 				 		tags$b("Choose the minimum and maximum spacer length: "),
						 				 		tags$br(),
						 				 		div(style="display:inline-block",
						 				 				numericInput("spamin",
						 				 										 label = "Min: ",
						 				 										 value = 14,
						 				 										 min = 14,
						 				 										 max = 16,
						 				 										 step = 1,
						 				 										 width = "80px"
						 				 				)),
						 				 		div(style="display:inline-block",
						 				 				numericInput("spamax",
						 				 										 label = "Max: ",
						 				 										 value = 16,
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
						 				 		
						 				 		textAreaInput("geneSeq",
						 				 									label = "",
						 				 									value = "",
						 				 									placeholder = "Paste gene sequence here..."),
						 				 		
						 				 		textOutput("validgeneseq"),
						 				 		p(""),
						 				 		
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
						 				 			rHandsontableOutput("exonInfo"),
						 				 			textOutput("validexoninfo")
						 				 		)
						 				 	)
						 				 ),
						 				 
						 				 
						 				 ############Exon Options#########################
						 				 conditionalPanel(
						 				 	condition = "input.inputType == 1",
						 				 	wellPanel(
						 				 		radioButtons("firstExon",
						 				 								 label = "2a. Do you want to find targets in the first exon? (Not recommended for gene knockouts.)",
						 				 								 choices = list("No" = 0,
						 				 								 							 "Yes" = 1),
						 				 								 selected = 0,
						 				 								 inline = TRUE),
						 				 		
						 				 		p(tags$b("Which exon(s) do you want to target?")),
						 				 		radioButtons("exonTargetType",
						 				 								 label = "Do you want to: ",
						 				 								 choices = list("Search all exons" = 0,
						 				 								 							 "Search for a target within a specified percentage of exons of the beginning of the sequence" = 1,
						 				 								 							 "Search for a target within a specified percentage of exons of the end of the sequence" = 2,
						 				 								 							 "Provide a list of exons to target" = 3),
						 				 								 selected = 0),
						 				 		
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
						 				 		
						 				 		conditionalPanel(
						 				 			condition = "input.exonTargetType == 2",
						 				 			p("Specify a percentage of exons, starting from the end of the sequence, to search. "),
						 				 			
						 				 			numericInput("exonEndPercentage",
						 				 									 label = "",
						 				 									 value = 30,
						 				 									 min = 1,
						 				 									 max = 100)
						 				 		),
						 				 		
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
						 				 	numericInput("threshold",
						 				 							 label = paste0("3. Choose minimum score for reporting (>50 is a strong score, ",
						 				 							 							 ">40 is moderate; we do not recommend using sites <40.): "),
						 				 							 value = 40,
						 				 							 min = 0,
						 				 							 max = 100
						 				 	),
						 				 	textOutput("validthreshold")
						 				 ),
						 				 
						 				 wellPanel(	
						 				 	
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 1",
						 				 		actionButton("genbankSubmit", "Submit"),
						 				 		p(""),
						 				 		
						 				 		uiOutput("downOutGenbank"),
						 				 		p(""),
						 				 		
						 				 		DT::dataTableOutput('genbankResults')
						 				 		#tags$head(tags$style("#geneSeqResults table {background-color: white; }", media = "screen", type = "text/css"))
						 				 	),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 2",
						 				 		actionButton("geneSeqSubmit", "Submit"),
						 				 		p(""),
						 				 		
						 				 		uiOutput("downOut"),
						 				 		p(""),
						 				 		
						 				 		DT::dataTableOutput('geneSeqResults')
						 				 		#tags$head(tags$style("#geneSeqResults table {background-color: white; }", media = "screen", type = "text/css"))
						 				 	)
						 				 )
						 				 )
						 ),
						 
						 
						 ##########HOW TO CITE Tab#########################################
						 tabPanel(
						 	tags$div("How to Cite", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		p("Paper link will be here when it's published.")
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		p("Manuscript is in prep; citation will be available shortly.")
						 	))
						 	
						 ),
						 
						 ##########CONTACT US Tab##########################################
						 tabPanel(
						 	tags$div("Contact", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		#p("Lab website links here.")
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		p("Please contact MENTHUHelp@gmail.com to report issues and request support."),
						 		p("Before submitting a bug report, please read the instructions below on how to write a helpful bug report."),
						 		p("By following these instructions, we will be able to solve your issue more quickly."),
						 		includeHTML("www/20170921_A_Guide_to_Writing_Helpful_Bug_Reports.html")
						 		
						 	))
						 	
						 ),
						 ##########TOOLS AND DOWNLOADS TAB#################################
						 
						 tabPanel(
						 	tags$div("Tools and Downloads", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		p(tags$a(href = "https://github.com/Dobbs-Lab/MENTHU", target = "_blank", "Download MENTHU at GitHub"))
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
						 		tags$html(tags$div(tags$span(a(href = "https://www.iastate.edu/", target = "_blank", "Iowa State University"))),
						 							tags$div(tags$span(a(href = "https://www.mayoclinic.org/", target = "_blank", "The Mayo Clinic"))),
						 							tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", target = "_blank", "Genome Writers Guild"))),
						 							tags$div(tags$span(a(href = "https://dill-picl.org", target = "_blank", "Dill-PICL Lab"))),
						 							tags$div(tags$span(a(href = "https://www.nih.gov/", target = "_blank", "National Institutes of Health (NIH)")))
						 		))),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		tags$p("This project was funded through NIH R24 OD020166 and is a joint effort by Iowa State University and The Mayo Clinic."),
						 		tags$p("This server is generously hosted by the", a(href = "https://dill-picl.org/", target = "_blank", 'Dill-PICL Lab'), "at Iowa State University.")
						 	))
						 ),
						 
						 
						 ######STATUS and CHANGELOG Tab####################################
						 tabPanel(
						 	tags$div("Change Log", style = "color:white"),
						 	titlePanel(""),
						 	
						 	column(2, wellPanel(
						 		tags$html(tags$div(tags$span(a(href = "https://www.iastate.edu/", target = "_blank", "Iowa State University"))),
						 							tags$div(tags$span(a(href = "https://www.mayoclinic.org/", target = "_blank", "The Mayo Clinic"))),
						 							tags$div(tags$span(a(href = "https://www.genomewritersguild.org/", target = "_blank", "Genome Writers Guild"))),
						 							tags$div(tags$span(a(href = "https://dill-picl.org", target = "_blank", "Dill-PICL Lab"))),
						 							tags$div(tags$span(a(href = "https://www.nih.gov/", target = "_blank", "National Institutes of Health (NIH)")))
						 		))),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 		includeHTML("www/changelog.html")
						 	))
						 )
	))
#