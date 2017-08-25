###############################################################################
#### * * * * * * * * * * * * * * * PIPELINE * * * * * * * * * * * * * * *  ####
###############################################################################
# This page contains detailed instructions on how to run the MENTHU pipeline.
# If you do not need detailed instructions, please see the "pipelineExpress.R"
# file for a pipeline with no instructions and only necessary fields.

####1. Replace the sequence below with your gene of interest####

yourGene <- "YOURGENESEQUENCEHERE"

####2. Use ONE of the options below:####
#  a. Choose a tab-delimited file containing the exon locations in your gene of 
#   interest; this file must be in the form:
#   EXON1START	EXON1END
#   EXON2START  EXON2END

exonArray <- read.table(file.choose(), stringsAsFactors = FALSE)

#  b. DO NOT USE THIS OPTION IF YOU USED OPTION A.
#   Replace the numbers below with the sequence indices of the exons in your
#   gene - e.g., if you have a gene with two exons located from nucleotide 189 to
#   322 and from 485 to 527, the code should look like this:
#   exonArray <- matrix(c(189, 322,
#											   485, 527),
#										   ncol = 2, 
#										   byrow = TRUE) 
#   Delete or add rows as necessary to list any exons you wish to target; it is
#   recommended that you choose an exon that is NOT the first exon, but codes for
#   the first 20-40% of the protein

exonArray <- matrix(c(1, 1,
											2, 2,
											3, 3,
											4, 4,
											5, 5),
										ncol = 2,     # don't change this
										byrow = TRUE) # don't change this

####3. OPTIONAL: Choose the tool type(s) that you will use to target the gene####
# If you do not change this, by default the code will search for listed tool
# type. You can put a '#' in front of options you don't want, e.g. if you do not
# want TALEN targets:
# targetList <- c("cas9" 
# 								#"talen",
# 								)

targetList <- c("cas9", 
								"talen"
								)


####4. Choose a name for the file that your results will be printed to####

outFile <- "YOURFILENAMEHERE.txt"


####5. Prediction step: Simply run this code after specifying your options; your####
# results will appear in the file name you specified in step 4.

menthu(yourGene, exonArray, targetList, outFile)




