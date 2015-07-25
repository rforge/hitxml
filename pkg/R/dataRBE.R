################################
# dataRBE CLASS
################################
setClass( Class            = "dataRBE",
          slots            = c( tissue           = "character",
                                alpha.1.Gy       = "numeric",
                                beta.1.Gy2       = "numeric",
                                D.cut.Gy         = "numeric",
                                r.nucleus.um     = "numeric",
                                RBE.alpha        = "data.frame"),
          prototype        = list(  tissue           = character(),
                                    alpha.1.Gy       = numeric(),
                                    beta.1.Gy2       = numeric(),
                                    D.cut.Gy         = numeric(),
                                    r.nucleus.um     = numeric(),
                                    RBE.alpha        = data.frame(E.MeV.u     = numeric(), 
                                                                  projectile  = character(),
                                                                  RBE.alpha   = numeric())))
################################
# Constructor
################################
dataRBE <- function(file.name, rbe.path = "."){
  
  # Nested function to extract header items
  read.item  <-  function(data, code){
    line  <-	input[grep(code, data)]
    if (!is.null(line)){
      return(as.numeric(substring( line, regexpr(" ", line) + 1, nchar(line))))
    }else{
      return(0.0)
    }
  }
  
  input		<-	scan(file.path(rbe.path,  
                           paste0(file.name, ".rbe"), 
                           fsep = .Platform$file.sep),
                   what = "character", 
                   sep = "\n")
  
  ######################
  # read projectile data
  lines			      <-	grep("!projectile", input)
  length.data.set	<-	diff(lines)[1] - 4		# Assumption - all datasets have same length, remove var names line
  df			        <-	NULL
  
  for(i in 1:length(lines)){
    #i <- 1
    projectile.name	<-	substring(input[lines[i]], 12, nchar(input[lines[i]]))
    projectile.name	<-	sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", projectile.name, perl=TRUE)
    xx			<-	input[(lines[i] + 3):(lines[i] + 3 + length.data.set)]
    xx			<-	sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", xx, perl=TRUE)
    tmp			<-	data.frame(	E.MeV.u 	= 	as.numeric(substring(xx, 1, regexpr(" ", xx) - 1)),
                            RBE.alpha =  	as.numeric(substring(xx, regexpr(" ", xx) + 1), nchar(xx)))
    tmp$projectile	<-	as.character(rep(projectile.name, nrow(tmp)))
    if(is.null(df)){
      df			<-	tmp
    }else{
      df			<-	rbind.data.frame(df, tmp)
    }
  }
  
  #############
  # read parameters
  alpha.1.Gy	  <-	read.item(input, "!alpha")
  beta.1.Gy2	  <-	read.item(input, "!beta")
  D.cut.Gy		  <-	read.item(input, "!cut")
  r.nucleus.um	<-	read.item(input, "!rnucleus")
  
  ############################
  # covert projectile to (A,Z)
  
  if(F){
    elements		<-	data.frame(	Z	= 1:10,
                           name	= c(	"H", "He", "Li", "Be", "B",
                                     "C", "N", "O", "F", "Ne"))
  
  df$element		<-	substring(df$projectile, regexpr("[[:alpha:]]", df$projectile, perl = TRUE), nchar(df$projectile))
  df$A			<-	as.numeric(substring(df$projectile, 1, regexpr("[[:alpha:]]", df$projectile, perl = TRUE) - 1))
  df$Z			<-	elements$Z[match(df$element, elements$name)]
  
  df$Z.A		<-	paste("(", sprintf("%02d", df$Z), ",", sprintf("%02d", df$A), ")", sep = "")
  }
  
  new("dataRBE",
      tissue       = gsub(".rbe", "", file.name),
      alpha.1.Gy   = alpha.1.Gy,
      beta.1.Gy2   = beta.1.Gy2,
      D.cut.Gy     = D.cut.Gy,
      r.nucleus.um = r.nucleus.um,
      RBE.alpha    = df)
}