################################
# dataDDD CLASS
################################
setClass( Class            = "dataDDD",
          slots            = c( projectile             = "character",
                                beam.energy.MeV.u      = "numeric",
                                target.material        = "character",
                                density.g.cm3          = "numeric",
                                peak.position.g.cm2    = "numeric",
                                DDD                    = "data.frame"),
          prototype        = list( projectile             = character(),
                                   beam.energy.MeV.u      = numeric(),
                                   target.material        = character(),
                                   density                = numeric(),
                                   peak.position.g.cm2    = numeric(),
                                   DDD                    = data.frame(depth.g.cm2                 = numeric(),
                                                                       dE.dz.MeV.cm2.g             = numeric())) )

################################
# Constructor
dataDDD <- function(file.name){

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

  ddd                         <- read.table( file.name, skip=10 )
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
# Method plot
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