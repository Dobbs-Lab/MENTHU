###############################################################################
#### * * * * * * * * * * * * * * * PIPELINE * * * * * * * * * * * * * * *  ####
###############################################################################
# If you need instructions for this pipeline, please see the file
# "pipelineInstructions.R"

yourGene <- "YOURGENESEQUENCEHERE"

exonArray <- read.table(file.choose(), stringsAsFactors = FALSE)

#exonArray <- matrix(c(1, 1,
#											2, 2,
#											3, 3,
#											4, 4,
#											5, 5),
#										ncol = 2,     
#										byrow = TRUE)

targetList <- c("cas9", 
								"talen"
)

outFile <- "YOURFILENAMEHERE.txt"

menthu(yourGene, exonArray, targetList, outFile)




