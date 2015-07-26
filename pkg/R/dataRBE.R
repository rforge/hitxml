################################
# dataRBE CLASS
################################
setClass( Class            = "dataRBE",
          slots            = c( cell.type        = "character",
                                file.date        = "character",
                                alpha.X          = "numeric",
                                beta.X           = "numeric",
                                D.cut.Gy         = "numeric",
                                r.nucleus.um     = "numeric",
                                RBE              = "data.frame"),
          prototype        = list(  cell.type        = character(),
                                    filedate         = character(),
                                    alpha.X          = numeric(),
                                    beta.X           = numeric(),
                                    D.cut.Gy         = numeric(),
                                    r.nucleus.um     = numeric(),
                                    RBE              = data.frame(E.MeV.u     = numeric(), 
                                                                  projectile  = character(),
                                                                  Z           = integer(),
                                                                  A           = integer(),
                                                                  RBE.initial = numeric())))
################################
# Constructor
################################
dataRBE <- function(file.name, rbe.path = "."){
  
  # Nested function to extract header items
  read.item.numeric  <-  function(data, code){
    line  <-	input[grep(code, data)]
    if (!is.null(line)){
      return(as.numeric(substring( line, regexpr(" ", line) + 1, nchar(line))))
    }else{
      return(0.0)
    }
  }
  read.item.character  <-  function(data, code){
    # returns string w/o leading or trailing whitespace
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    
    line  <-	input[grep(code, data)]
    if (!is.null(line)){
      return(trim(substring( line, regexpr(" ", line) + 1, nchar(line))))
    }else{
      return(0.0)
    }
  }

  input		<-	scan(file.path(rbe.path,  
                             file.name, 
                             fsep = .Platform$file.sep),
                   what = "character", 
                   sep = "\n")
  
  #################
  # read parameters
  file.type     <-  read.item.character(input, "!filetype")
  if(file.type != "RBE"){
    stop("File is not a valid RBE data file.")
  }
  file.date     <-  read.item.character(input, "!filedate")
  alpha.X    	  <-	read.item.numeric(input, "!alpha")
  beta.X    	  <-	read.item.numeric(input, "!beta")
  D.cut.Gy		  <-	read.item.numeric(input, "!cut")
  r.nucleus.um	<-	read.item.numeric(input, "!rnucleus")
  
  
  ######################
  # read projectile data
  lines			      <-	grep("!projectile", input)
  length.data.set	<-	diff(lines)[1] - 4		# Assumption - all datasets have same length, remove var names line
  df			        <-	NULL
  
  for(i in 1:length(lines)){
    #i <- 1
    projectile.name	<-	substring(input[lines[i]], 12, nchar(input[lines[i]]))
    projectile.name	<-	sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", projectile.name, perl=TRUE)
    xx			        <-	input[(lines[i] + 3):(lines[i] + 3 + length.data.set)]
    xx			        <-	sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", xx, perl=TRUE)
    tmp			        <-	data.frame(	E.MeV.u 	  = 	as.numeric(substring(xx, 1, regexpr(" ", xx) - 1)),
                                    RBE.initial =  	as.numeric(substring(xx, regexpr(" ", xx) + 1), nchar(xx)))
    tmp$projectile	<-	as.character(rep(projectile.name, nrow(tmp)))
    particle.no     <-  AT.particle.no.from.particle.name(tmp$projectile)
    tmp$Z           <-  AT.Z.from.particle.no(particle.no)$Z
    tmp$A           <-  AT.A.from.particle.no(particle.no)$A
    if(is.null(df)){
      df			<-	tmp
    }else{
      df			<-	rbind.data.frame(df, tmp)
    }
  }
  

  new("dataRBE",
      cell.type    = gsub(".rbe", "", file.name),
      file.date    = file.date,
      alpha.X      = alpha.X,
      beta.X       = beta.X,
      D.cut.Gy     = D.cut.Gy,
      r.nucleus.um = r.nucleus.um,
      RBE          = df)
}


################################
# Getter functions
################################
alpha.X <- function(x){
  if(class(x) != "dataRBE"){
    stop("x must be of class 'dataRBE'")
  }
  if(is.null(x@alpha.X)){
    stop("No value for alpha found")
  }else{
    return(x@alpha.X)
  }
}

beta.X <- function(x){
  if(class(x) != "dataRBE"){
    stop("x must be of class 'dataRBE'")
  }
  if(is.null(x@beta.X)){
    stop("No value for beta found")
  }else{
    return(x@beta.X)
  }
}

D.cut.Gy <- function(x){
  if(class(x) != "dataRBE"){
    stop("x must be of class 'dataRBE'")
  }
  if(is.null(x@D.cut.Gy)){
    stop("No value for Dcut found")
  }else{
    return(x@D.cut.Gy)
  }
}

r.nucleus.um <- function(x){
  if(class(x) != "dataRBE"){
    stop("x must be of class 'dataRBE'")
  }
  if(is.null(x@r.nucleus.um)){
    stop("No value for nucleus radius found")
  }else{
    return(x@r.nucleus.um)
  }
}

#  Computes the maximum slope of the X-ray dose-effect curve
s.max  <- function(x){
  return(alpha.X(x) + (2 * beta.X(x) * D.cut.Gy(x)) )
}


# FUNCTION RBE.initial
#  Returns initial RBE (values stored in RBE file) for
#  a series of projektiles and energies
#
# Arguments
#  x            dataRBE object
#  projectile   particle name
#  E.MeV.u      energy of the particle
#
RBE.initial <- function(x, projectile, E.MeV.u){
  if(length(projectile) != length(E.MeV.u)){
    stop("Arguments must have same length.")
  }
  
  RBE   <- numeric(length(projectile))
  
  for(cur.projectile in unique(projectile)){
    # cur.projectile <- unique(projectile)[2]
    ii        <- projectile == cur.projectile
    jj        <- x@RBE$projectile == cur.projectile
    no.data   <- FALSE
    if(sum(jj) == 0){
      Z <- AT.Z.from.particle.no(AT.particle.no.from.particle.name(cur.projectile))$Z
      fall.back.projectile <- AT.particle.name.from.particle.no(Z*1000+2*Z)
      warning(paste0("Projectile ", 
                     cur.projectile,
                     " not available in RBE data. Trying ",
                     fall.back.projectile,
                     " instead."))
      jj <- x@RBE$projectile == fall.back.projectile
      if(sum(jj) == 0){
        warning(paste0("Fall back projectile ", 
                       cur.projectile,
                       " also not available in RBE data. Returning NAs."))
        no.data <- TRUE
      }
    }
    if(no.data){
      RBE[ii] <- NA
    }else{
      RBE[ii]   <- approx( x@RBE$E.MeV.u[jj], 
                           x@RBE$RBE.initial[jj], 
                           xout=E.MeV.u[ii])$y
    }
  }
  return(RBE)
}


# Computes the initial slope of the ion dose-effect curve
#
# Arguments
#  x            dataRBE object
#  projectile   particle name
#  E.MeV.u      energy of the particle
#
alpha.ion  <- function(x, projectile, E.MeV.u){
  return(alpha.X(x) * RBE.initial(x, projectile, E.MeV.u))
}

# Computes the beta coefficient for the ion dose-effect curve
#
# Arguments
#  x            dataRBE object
#  projectile   particle name
#  E.MeV.u      energy of the particle
#
beta.ion  <- function(x, projectile, E.MeV.u){
  return( s.max(x) - alpha.ion(x, projectile, E.MeV.u) / (2 * D.cut.Gy(x)) )
}

