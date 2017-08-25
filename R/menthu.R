
menthu <- function(seq, exon, targetList, outFile){
	#Identify all potential targets
	targets <- getTargets(seq, exon, targetList)
	
	#Calculate Bae and competition score at all potential targets
	scores <- calculateScores(targets)
	
	#Generate targeting sequence for targets?
	
	#Sort and format outputs
	outFrame <- data.frame(target_gene_id         = targets$id,
												 target_sequence        = targets$seq,
												 Bae_out-of-frame_score = scores$oofScore,
												 menthuScore            = scores$menthuScore,
												 tool                   = targets$toolType,
												 guide_sequence         = ,
												 stringsAsFactors = FALSE
												 )
	
	
	
	
	#Output outputs
	printMenthu(menthuScores)
}