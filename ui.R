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
#library(rhandsontable)
#library(shinyIncubator)

shinyUI(
	
	
	#Creates the navbar set-up
	navbarPage(id = 'mainPage',
						 windowTitle = "MENTHU",
						 
						 #Stylesheet
						 theme = "ogtheme.css", 
						 
						 #Page title box
						 tags$div("MENTHU v1.0.0", 
						 				 style = "color:white"),
						 
						 ########ABOUT TAB#################################################
						 tabPanel(tags$div("About", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isulogo.jpg",     width = "90%")))),
						 				 	p(""),
						 				 	tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "mayoclinic.jpeg", width = "50%"))))
						 				 )),
						 				 
						 				 #Text area in center of page
						 				 column(9, wellPanel(
						 				 ))
						 				 
						 ),
						 
						 ##########INSTRUCTIONS############################################
						 tabPanel(tags$div("Instructions and FAQs", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 #Sidebar panel with links
						 				 column(2, wellPanel(
						 				 	tags$div(tags$span(a(href = "https://www.iastate.edu/",   target = "_blank", tags$img(src = "isulogo.jpg",     width = "90%")))),
						 				 	p(""),
						 				 	tags$div(tags$span(a(href = "https://www.mayoclinic.org", target = "_blank", tags$img(src = "mayoclinic.jpeg", width = "50%"))))
						 				 )),
						 				 
						 				 #Text area in center of page
						 				 column(9, wellPanel(
						 				 	
						 				 	#Display page instructions
						 				 	#includeHTML("www/instructions.html")
						 				 ))
						 				 
						 ),
						 
						 ##########SUBMIT SINGLE JOB TAB###################################
						 tabPanel(id = "single",
						 				 tags$div("Calculate MENTHU Score", style = "color:white"),
						 				 titlePanel(""),
						 				 
						 				 ##Sidebar############################################################
						 				 #Adds a sidebar for users to pre-populate fields with an example, and reset the form
						 				 column(2, wellPanel(
						 				 	
						 				 	#Cut/Paste cDNA example; input$example
						 				 	actionLink("example", 
						 				 						 label = "Example"),
						 				 	
						 				 	p(""),
						 				 	
						 				 	#Reset Button; input$reset
						 				 	actionLink("reset", 
						 				 						 label = "Reset Form")
						 				 )),
						 				 
						 				 
						 				 ####Main Bar#########################################################
						 				 #Main panel for entering information and submitting job
						 				 column(9, wellPanel(
						 				 	
						 				 	####Choose PAM sequence#############################################
						 				 	p("1. Select the PAM sequence(s) you wish to target:"),
						 				 	
						 				 	actionLink("selectAll", "Select All"),
						 				 	p(""),
						 				 	actionLink("selectNone", "De-select All"),
						 				 	checkboxGroupInput("casType",
						 				 										 label = "",
						 				 										 choices = list("S. pyogenes SpCas9: 5'-NGG-3'" = 1,
						 				 										 							 "S. pyogenes SpCas9: 5'-NRG-3'" = 2,
						 				 										 							 "S. aureus SaCas9: 5'-NNNRRT-3'" = 3,
						 				 										 							 "S. aureus SaCas9: 5'-NNGRRT-3'" = 4,
						 				 										 							 "S. pasteurianus SpCas9: 5'-NNGTGA-3'" = 5,
						 				 										 							 "S. thermophilus StCas9: 5'-NNAGAAW-3'" = 6,
						 				 										 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTN-3'" = 7,
						 				 										 							 #"Acidaminococcus AsCpf1/Lachnospiraceae LbCpf1: 5'-TTTV-3'" = 8,
						 				 										 							 "C. jejuni CjCas9: 5'-NNNVRYAC-3'" = 9,
						 				 										 							 #"Francisella FnCpf1: 5'-TTN-3'" = 10,
						 				 										 							 #"Francisella FnCpf1: 5'-YTN-3'" = 11,
						 				 										 							 "N. meningitidis NmCas9: 5'-NNNNGMTT-3'" = 12),
						 				 										 							 #"T. denticola: 5'-NAAAC-3'" = 13),
						 				 										 selected = 1)
						 				 ),
						 				 wellPanel(
						 				 	
						 				 	
						 				 	
						 				 	####Choose Input Type###############################################
						 				 	p("2. What type of sequence input do you want to use?"),
						 				 	radioButtons("inputType",
						 				 							 label = "",
						 				 							 choices = list(#"GenBank Gene ID"  = 1,
						 				 							 							 "Copy/paste sequence" = 2
						 				 							 ),
						 				 							 selected = 1),
						 				 	####Input sequence#####################
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 1",
						 				 		
						 				 		textAreaInput("genbankId",
						 				 									label = "",
						 				 									value = "",
						 				 									placeholder = "Paste GenBank gene ID here...")
						 				 		
						 				 		
						 				 	),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 2",
						 				 		
						 				 		textAreaInput("geneSeq",
						 				 									label = "",
						 				 									value = "",
						 				 									placeholder = "Paste gene sequence here..."),
						 				 		p(""),
						 				 		p("Input exon location information here."),
						 				 	  p("If you need to add more exons, right-click and select 'Insert Row'")#,
						 				 		#rHandsontableOutput("exonInfo")
						 				 		
						 				 	)
						 				 	
						 				 ),
						 				 
						 				 
						 				 ############Exon Options#########################
						 				 
						 				 wellPanel(
						 				 	p("3. Do you want to use targets in the first exon? (Not recommended.)"),
						 				 	radioButtons("firstExon",
						 				 							 label = "",
						 				 							 choices = list("No" = 0,
						 				 							 							 "Yes" = 1),
						 				 							 selected = 0)
						 				 ),
						 				 
						 				 wellPanel(
						 				 	p("4. What percentage of the first exons in the sequence do you want to consider?"),
						 				 	p("For gene knockouts, it is recommended that you find a target within the first 30% of the exons in the sequence."),
						 				 	p("However, you can look for targets in every exon (100%) by changing the value to '100'."),
						 				 	
						 				 	numericInput("exonPercentage",
						 				 							 label = "",
						 				 							 value = "30")
						 				 	
						 				 	
						 				 ),
						 				 
						 				 wellPanel(	
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 1",
						 				 		actionButton("genBankSubmit", "Submit")
						 				 	),
						 				 	
						 				 	conditionalPanel(
						 				 		condition = "input.inputType == 2",
						 				 		actionButton("geneSeqSubmit", "Submit")
						 				 	)
						 				 	
						 				 	
						 				 	
						 				 ))
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
						 	column(9,
						 				 
						 				 wellPanel(
						 		p("Please contact MENTHUHelp@gmail.com to report issues and request support."),
						 		p("Before submitting a bug report, please read the instructions below on how to write a helpful bug report."),
						 		p("By following these instructions, we will be able to solve your issue more quickly.")),
						 		wellPanel(
						 		includeHTML("www/20170921_A_Guide_to_Writing_Helpful_Bug_Reports.html"))
						 		)
						 	
						 ),
						 ##########TOOLS AND DOWNLOADS TAB#################################
						 
						 tabPanel(
						 	tags$div("Tools and Downloads", style = "color:white"),
						 	titlePanel(""),
						 	#Sidebar panel with links
						 	column(2, wellPanel(
						 		p(tags$a(href = "https://github.com/Dobbs-Lab/GTagHD", target = "_blank", "Download GTagHD at GitHub"))
						 	)),
						 	
						 	#Text area in center of page
						 	column(9, wellPanel(
						 			p("A standalone version of this code will be available shortly at GitHub.")
						 			#p("A standalone version of this code may be downloaded from", tags$a(href = "https://github.com/Dobbs-Lab/GTagHD", target = "_blank", " GitHub."), " The R code is provided as-is, and may not be used in commercial applications. Please be aware that you modify the code at your own risk; we are unable to provide support for modified versions.")
						 		
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
						 		#includeHTML("www/changelog.html")
						 	))
						 )
	))
