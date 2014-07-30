leading.blanks <- function(df, blank.char = ".", extra.char = " ", min.width = 1, max.width = 99, ignore.first = F, trace = F){
	# Created: April 10, 2005
	# Revised: April 10, 2005
	# Name   : Claus E. Andersen
	# Main purpose: to align data in Latex tables.
	# The ignore.first parameters controls the first row (often a row with names).
	# For example:
	#   xx <- data.frame(sample(c(1,1055500,10),10,replace=T),sample(c("a","aba","sfgsdfgfhgd"),10,replace=T))
	#   xx.LATEX <- create.latex.table(xx)
	#   yy <- leading.blanks(xx.LATEX,"x")
	#   yy
	#   write.table(yy,"xx.txt",sep=" ",dimnames.write = F)
	for(i in 1:ncol(df)) {
		cc <- df[, i]
		cc <- as.character(cc)
		cc.all <- cc
		if(ignore.first) cc <- cc[-1]
		cc.width <- nchar(as.character(cc)) #S-plus: string.bounding.box(cc)$columns
		cc.max.width <- min(max.width, max(cc.width, min.width))
		cc.pure <- paste(rep(blank.char, cc.max.width), sep = "",
			collapse = "")
		if(trace)
			print(paste("Col = ", i, " class = ", class(cc), 
				cc.max.width))
		if(ignore.first) {
			df[, i] <- c(cc.all[1], paste(extra.char, paste(
				substring(cc.pure, 1, cc.max.width - cc.width),
				cc, sep = ""), sep = ""))
		}
		else {
			df[, i] <- paste(extra.char, substring(cc.pure, 1,
				cc.max.width - cc.width), cc, sep = "")
		}
	}
	df
}

#########################################################
 
