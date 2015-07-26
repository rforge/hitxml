HX.read.RBE.file <- function(file.name){
  
  # Nested function to extract header items
  read.item  <-  function(data, code){
    line  <-	input[grep(code, data)]
    if (!is.null(line)){
      return(as.numeric(substring( line, regexpr(" ", line) + 1, nchar(line))))
    }else{
      return(0.0)
    }
  }
  
  # DEBUG file.name <- file.names[1]
  input		<-	scan(file.name, what = "character", sep = "\n")
  
  
  ######################
  # read projectile data
  lines			<-	grep("!projectile", input)
  length.data.set	<-	diff(lines)[1] - 4		# Assumption - all datasets have same length, remove var names line
  df			<-	NULL
  
  for(i in 1:length(lines)){
    #i <- 1
    projectile.name	<-	substring(input[lines[i]], 12, nchar(input[lines[i]]))
    projectile.name	<-	sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", projectile.name, perl=TRUE)
    xx			<-	input[(lines[i] + 3):(lines[i] + 3 + length.data.set)]
    xx			<-	sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", xx, perl=TRUE)
    tmp			<-	data.frame(	E.MeV.u 	= 	as.numeric(substring(xx, 1, regexpr(" ", xx) - 1)),
                         RBE		=  	as.numeric(substring(xx, regexpr(" ", xx) + 1), nchar(xx)))
    tmp$projectile	<-	as.character(rep(projectile.name, nrow(tmp)))
    if(is.null(df)){
      df			<-	tmp
    }else{
      df			<-	rbind.data.frame(df, tmp)
    }
  }
  
  #############
  # add header
  df$alpha.1.Gy	<-	rep(read.item(input, "!alpha"), nrow(df))
  df$beta.1.Gy2	<-	rep(read.item(input, "!beta"), nrow(df))
  df$D.cut.Gy		<-	rep(read.item(input, "!cut"), nrow(df))
  df$r.nucleus.um	<-	rep(read.item(input, "!rnucleus"), nrow(df))
  
  ############################
  # covert projectile to (A,Z)
  
  elements		<-	data.frame(	Z	= 1:10,
                           name	= c(	"H", "He", "Li", "Be", "B",
                                     "C", "N", "O", "F", "Ne"))
  
  df$element		<-	substring(df$projectile, regexpr("[[:alpha:]]", df$projectile, perl = TRUE), nchar(df$projectile))
  df$A			<-	as.numeric(substring(df$projectile, 1, regexpr("[[:alpha:]]", df$projectile, perl = TRUE) - 1))
  df$Z			<-	elements$Z[match(df$element, elements$name)]
  
  df$Z.A		<-	paste("(", sprintf("%02d", df$Z), ",", sprintf("%02d", df$A), ")", sep = "")
  
  return(df)
}
