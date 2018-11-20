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
                                   alpha.X.Gy             = numeric(),
                                   beta.X.Gy2             = numeric(),
                                   RBE.model              = character(),
                                   DDD                    = data.frame(depth.g.cm2                 = numeric(),
                                                                       dE.dz.MeV.cm2.g             = numeric(),
                                                                       alpha.ion.Gy                = numeric(),
                                                                       beta.ion.Gy2                = numeric())) )

################################
# Constructor
dataDDD <- function(file.name){
  # file.name <- x

  # Nested function to extract header items
  read.item.numeric  <-  function(data, code){
    line  <-	input[grep(code, data)]
    if (length(line) != 0){
      line <- gsub("\\s+", " ", stringr::str_trim(line))
      parts <- strsplit(line, " ")
      return(as.numeric(parts[[1]][2]))
    }else{
      return(NA)
    }
  }
  
  read.item.character  <-  function(data, code){
    # returns string w/o leading or trailing whitespace
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    
    line  <-	input[grep(code, data)]
    if (!is.null(line)){
      return(trim(substring( line, regexpr(" ", line) + 1, nchar(line))))
    }else{
      return("N/A")
    }
  }
  
  # Read in TRiP98/syngoRT ddd files, potential extended by rbe data
  input		<-	scan(file.name, 
                 what = "character", 
                 sep = "\n")
  
  #################
  # read parameters
  file.type     <-  read.item.character(input, "!filetype")
  if(!file.type%in%c("ddd", "ddd.HITXML")){
    stop("File is not a valid DDD data file.")
  }
  file.date     <-  read.item.character(input, "!filedate")
  projectile	  <-	read.item.character(input, "!projectile")
  material    	<-	read.item.character(input, "!material")
  density  		  <-	read.item.numeric(input, "!density")
  energy      	<-	read.item.numeric(input, "!energy")
  alpha.X.Gy    <-	read.item.numeric(input, "!alphaX:")
  beta.X.Gy2    <-	read.item.numeric(input, "!betaX:")
  RBE.model     <-  read.item.character(input, "!rbemodel")
  start.line    <-  grep("!ddd", input)[1]
  ddd           <-  read.table( skip = start.line,
                                sep = " ",
                                text = gsub(",", " ",
                                       gsub("\t", " ", 
                                            readLines(file.name))))

  if(file.type == "ddd"){
    if(ncol(ddd)>2){
      ddd <- ddd[,seq(-3, -1*ncol(ddd))]
    }
    ddd[3] <- rep(NA, nrow(ddd))
    ddd[4] <- rep(NA, nrow(ddd))
  }
  
  colnames(ddd)               <- c("depth.g.cm2", "dE.dz.MeV.cm2.g", "alpha.ion.Gy", "beta.ion.Gy2")

  new("dataDDD",
      projectile          = projectile,
      beam.energy.MeV.u   = energy,
      target.material     = material,
      density.g.cm3       = density,
      peak.position.g.cm2 = ddd$depth.g.cm2[which.max(ddd$dE.dz.MeV.cm2.g)],
      alpha.X.Gy          = alpha.X.Gy,
      beta.X.Gy2          = beta.X.Gy2,
      RBE.model           = RBE.model,
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

get.alpha.ion.Gy <- function(DDD.data, depths.g.cm2){
  alphas <- approx(DDD.data@DDD$depth.g.cm2, 
                  DDD.data@DDD$alpha.ion.Gy, 
                  xout = depths.g.cm2)$y
  alphas[is.na(alphas)] <- 0.0
  
  return( alphas )
}

get.beta.ion.Gy2 <- function(DDD.data, depths.g.cm2){
  betas <- approx(DDD.data@DDD$depth.g.cm2, 
                   DDD.data@DDD$beta.ion.Gy2, 
                   xout = depths.g.cm2)$y
  betas[is.na(betas)] <- 0.0
  
  return( betas )
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
