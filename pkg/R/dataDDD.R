################################
# dataDDD CLASS
################################
setClass( Class            = "dataDDD",
          slots            = c( projectile             = "character",
                                beam.energy.MeV.u      = "numeric",
                                target.material        = "character",
                                density.g.cm3          = "numeric",
                                peak.position.g.cm2    = "numeric",
                                alpha.X.Gy             = "numeric",
                                beta.X.Gy2             = "numeric",
                                RBE.model              = "character",
                                DDD                    = "data.frame"),
          prototype        = list( projectile             = character(),
                                   beam.energy.MeV.u      = numeric(),
                                   target.material        = character(),
                                   density                = numeric(),
                                   peak.position.g.cm2    = numeric(),
                                   DDD                    = data.frame(depth.g.cm2                 = numeric(),
                                                                       dE.dz.MeV.cm2.g             = numeric(),
                                                                       alpha.ion.Gy                = numeric(),
                                                                       beta.ion.Gy2                = numeric())) )

################################
# Constructor
dataDDD <- function(file.name, type){
  # file.name <- x

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
  
  if(!tools::file_ext(file.name)%in%c("ddd", "csv")){
    stop("Error: unknown DDD file format")
  }
  
  # Read in TRiP98/syngoRT ddd files
  if(tools::file_ext(file.name) == "ddd"){
    input		<-	scan(file.name, 
                 what = "character", 
                 sep = "\n")
  
    #################
    # read parameters
    file.type     <-  read.item.character(input, "!filetype")
    if(file.type != "ddd"){
      stop("File is not a valid DDD data file.")
    }
    file.date     <-  read.item.character(input, "!filedate")
    projectile	  <-	read.item.character(input, "!projectile")
    material    	<-	read.item.character(input, "!material")
    density  		  <-	read.item.numeric(input, "!density")
    energy      	<-	read.item.numeric(input, "!energy")
    
    start.line    <- grep("!ddd", input)[1]
    ddd                         <- read.table( file.name, skip=start.line )
  }
  
  # Read in Lucas style csv files
  if(tools::file_ext(file.name) == "csv"){
    input		<-	scan(file.name, 
                   what = "character", 
                   sep = "\n")
    
    #################
    # read parameters
    file.type     <-  "csv"
    file.date     <-  "n/a"
    projectile	  <-	"n/a"
    material    	<-	"n/a"
    density  		  <-	"n/a"
    energy      	<-	read.item.numeric(input, "# nominal energy: ")
    
    ddd                         <- read.table( file.name, skip=10 )
  }
  
  
  
  
  if(ncol(ddd)>2) {
      ddd <- ddd[,-5]
      ddd <- ddd[,-4]
      ddd <- ddd[,-3]
    }
  
  colnames(ddd)               <- c("depth.g.cm2", "dE.dz.MeV.cm2.g")

  new("dataDDD",
      projectile          = projectile,
      beam.energy.MeV.u   = energy,
      target.material     = material,
      density.g.cm3       = density,
      peak.position.g.cm2 = ddd$depth.g.cm2[which.max(ddd$dE.dz.MeV.cm2.g)],
      DDD                 = ddd)
}

get.dose.Gy <- function(DDD.data, depths.g.cm2){
  doses <- approx(DDD.data@DDD$depth.g.cm2, 
                  DDD.data@DDD$dE.dz.MeV.cm2.g, 
                  xout = depths.g.cm2)$y *
           DDD.data@density.g.cm3 * 1.6022e-10
  doses[is.na(doses)] <- 0.0
  
  return( doses )
}

######################
# Methods
setMethod(f          = "plot", 
          signature  = c("dataDDD"),
          definition = function(x) {
            
            lattice::xyplot(dE.dz.MeV.cm2.g ~ depth.g.cm2,
                            x@DDD,
                            type     = "o",
                            grid     = TRUE,
                            ylab     = "dE/dz / (MeV*cm2/g)",
                            xlab     = "depth / (g / cm2)",
                            main     = paste0(x@projectile,
                                              " (",
                                              x@beam.energy.MeV.u,
                                              " MeV/u) on ",
                                              x@target.material,
                                              " - peak at ",
                                              format(x@peak.position.g.cm2, digits = 3),
                                              " g/cm2"))
          })

setMethod(f          = "*", 
          signature  = c("dataDDD", "numeric"),
          definition = function(e1, e2) {
          
            new.ddd                  <- e1@DDD
            new.ddd$dE.dz.MeV.cm2.g  <- new.ddd$dE.dz.MeV.cm2.g * e2
            
            new("dataDDD",
                projectile          = e1@projectile,
                beam.energy.MeV.u   = e1@beam.energy.MeV.u,
                target.material     = e1@target.material,
                density.g.cm3       = e1@density.g.cm3,
                peak.position.g.cm2 = e1@peak.position.g.cm2,
                DDD                 = new.ddd)
          })

######################
# Routines
writeDDD <- function(ddd){
  output <- c("!filetype    ddd", 
              "!fileversion    19980520", 
              paste0("!filedate    ", date()), 
              paste0("!projectile    ", ddd@projectile),
              "!material      H2O",
              "!composition   H2O",
              "!density 1",
              paste0("!energy ", ddd@beam.energy.MeV.u),
              "   z[g/cm**2] dE/dz[MeV/(g/cm**2)] FWHM1[g/cm**2] factor FWHM2[g/cm**2]\n!ddd\n")
  write(output, "test.ddd", sep = "\n")
}
