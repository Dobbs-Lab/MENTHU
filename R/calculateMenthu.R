


calculateMenthu <- function(baeFrame){
	#Subset the Bae patterns to only include those where microhomology arm length is 3 or greater
	menthuPats <- baeFrame[which(nchar(baeFrame$microhomology) > 2),]
	#Order the patterns by pattern score
	menthuPats <- menthuPats[order(-menthuPats$patternScore),]
	
	#Create a data frame containing the top 10 or less pattern scores
	if(nrow(menthuPats) < 10){
		end <- nrow(menthuPats)
	} else {
		end <- 10
	}
	
	#Subset data to top 10
	top10 <- menthuPats[1:end,]
	
	#Perform linear regression on the top10 data
	linModel <- lm(top10$patternScore ~ seq(1:end), data = top10)
	
	return(linModel$coefficients[2])
}

