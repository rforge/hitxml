spectra <- function(nrow){
  return(matrix(nrow     = nrow,
                ncol     = 6,
                dimnames = list(NULL,
                                c("depth.step",
                                  "depth.g.cm2",
                                  "particle.no",
                                  "E.MeV.u",
                                  "dE.MeV.u",
                                  "N.per.primary"))))
}

################################
# dataSPC CLASS
setClass( Class            = "dataSPC",
          slots            = c( projectile             = "character",
                                beam.energy.MeV.u      = "numeric",
                                target.material        = "character",
                                peak.position.g.cm2    = "numeric",
                                redistributed          = "logical",
                                spectra                = "matrix"),
          prototype        = list( projectile             = character(),
                                   beam.energy.MeV.u      = numeric(),
                                   target.material        = character(),
                                   peak.position.g.cm2    = numeric(),
                                   redistributed          = logical(),
                                   spectra                = spectra(0)))


################################
# Constructor
dataSPC <- function(file.name, 
                    endian       = c("big", "little")[1], 
                    redistribute = TRUE,
                    ...){
  res   <- SPC.read( file.name = file.name,
                     endian    = endian)
  # Replace target material name by libamtrack compatible
  if(res$target.material == "H2O"){
    res$target.material = "Water, Liquid"
  }

  spc           <- as.matrix(res$spc[,!(names(res$spc)%in%c("E.low.MeV.u", 
                                                            "E.high.MeV.u", 
                                                            "dN.dE.per.MeV.u.per.primary"))])

  new.spc <- new( "dataSPC",
                  projectile          = res$projectile,
                  beam.energy.MeV.u   = res$energy.MeV.u,
                  target.material     = res$target.material,
                  peak.position.g.cm2 = res$peak.position.g.cm2,
                  redistributed       = FALSE, 
                  spectra             = spc)
  if(redistribute){
    new.spc <- redistribute.spc(new.spc, ...)}
  
  return(new.spc)
  }

# dataSPC <- function(projectile, beam.energy.MeV.u, target.material,
#                     peak.position.g.cm2, depth.step, depth.g.cm2,
#                     particle.no, E.MeV.u, dE.MeV.u, N.per.primary){
#   return( new( "dataSPC",
#                projectile          = projectile,
#                beam.energy.MeV.u   = beam.energy.MeV.u,
#                target.material     = target.material,
#                peak.position.g.cm2 = peak.position.g.cm2,
#                redistributed       = FALSE, 
#                spectra             = matrix(data     = c(depth.step,
#                                                          depth.g.cm2,
#                                                          particle.no,
#                                                          E.MeV.u,
#                                                          dE.MeV.u,
#                                                          N.per.primary),
#                                             nrow     = length(depth.step),
#                                             ncol     = 6,
#                                             dimnames = list(NULL,
#                                                             c("depth.step",
#                                                               "depth.g.cm2",
#                                                               "particle.no",
#                                                               "E.MeV.u",
#                                                               "dE.MeV.u",
#                                                               "N.per.primary")))))
# }

################################
# Get depth-dose
depth.dose <- function(x){
  DDD           <- SPC.tapply( x                    = x, 
                               INDEX                = c("depth.g.cm2"), 
                               FUN                  = AT.total.D.Gy, 
                               additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                           c("stopping.power.source.no", "2", FALSE)),
                               names.results        = "D.Gy")

  # Translate to ddd
  conversion.factor <- 6.2415e9
  DDD$D.Gy          <- DDD$D.Gy * conversion.factor
  return(DDD)
}

#####################################
# R function for reading of spc files 
#####################################
#
# started 2010/02/18, sgre
# as R script
#
# revised 2010/05/18, fk
# added: show.bin/show.bin.width feature, set type of mean.
# Bug fixed: tmp.data$Cum changed to tmp.data$HCum after tmp.data$H.
#
# revised 2010/07/28, sgre
# removed DDD handling
#
# revised 2010/08/04, sgre
# now choice of endian
#
# revised 2010/11/02, sgre
# added to R package
#
# revised 2010/11/10, sgre
# switched to matrix ops instead of dataframes, optimized for libamtrack, removed raw option
#
# revised 2010/11/11, sgre
# speed-up due to working on index tables, prototype for C version
#
# revised 2011/03/12, sgre
# after being superseeded by C version (r781) which
# still at troubles at r904 rolled back and combined
# with C version
# Now both R ('vanilla') and C (unstable) version can be used
#
# revised 2011/09/19, sgre
# removed "compress option" since it interferes with
# interpolation (maybe not in R but in later C translation of this code)
#
# revised 2015/07/26, sgre
# Moved from libamtrack to HITXML
#
SPC.read <- function( file.name, 
                      flavour        = "vanilla",
                      endian         = c("big", "little")[1], 
                      mean           = c("geometric", "arithmetic")[2],
                      header.only    = FALSE)
{
  # Unzip file, if applicable
  file.name.org <- file.name
  file.name     <- sub("\\.[[:alnum:]]+$", "", basename(as.character(file.name)))
  file.ext      <- substring(file.name.org, nchar(file.name)+2)
  if(file.ext == "zip"){
    unzip(file.name.org)
    file.name   <- paste(file.name, ".spc", sep = "")
  }else{
    file.name   <- file.name.org
  }
  
  if(flavour == "vanilla"){
    # R version
    to.read                 <- file(file.name, "rb")
    
    # Scan file first to get bin sizes
    seek(to.read, where = 0, origin = 'end')  
    end.pos                 <- seek(to.read, where = NA)  
    seek(to.read, where = 0, origin = 'start')
    tags                    <- 1
    max.tags                <- 20000
    mtag                    <- matrix( data = 0, ncol = 4, nrow = max.tags)
    repeat{
      mtag[tags, 1]  <- seek(to.read, where = NA)  
      mtag[tags, 2]  <- readBin(to.read, integer(), endian = endian)
      if(mtag[tags, 2] > 20){
        cat("Strange data read. Probably wrong endianess given.\n")
        close(to.read)
        stop()
      }
      mtag[tags, 3]  <- readBin(to.read, integer(), endian = endian)
      if( mtag[tags, 2] %in% c(9,12,16,18)){
        mtag[tags, 4]  <- readBin(    to.read, integer(), size = 8, n = 1, endian = endian)
      }else if(mtag[tags, 2] == 10){
        mtag[tags, 4]  <- readBin(    to.read, double(), size = 8, n = 1, endian = endian)
      }else if(mtag[tags, 2] == 13){
        seek(to.read, where = 16, origin = 'current')
        mtag[tags, 4]  <- 1000 * readBin( to.read, integer(), size = 4, signed = TRUE, n = 1, endian = endian) + readBin( to.read, integer(), size = 4, signed = TRUE, n = 1, endian = endian)
      } else{
        seek(to.read, where = mtag[tags, 3], origin = 'current')
      }
      tags <- tags + 1
      if(tags >= max.tags){
        cat("Exceeded maximum number of tags. Please choose higher number.\n")
        close(to.read)
        stop()
      }
      if(seek(to.read, where = NA) >= end.pos){
        break()
      }
    }
    
    mtag               <- mtag[mtag[,1] != 0,]
    n.depth.steps      <- mtag[mtag[,2] == 9, 4]
    depth.g.cm2        <- mtag[mtag[,2] == 10, 4]
    n.particle.species <- mtag[mtag[,2] == 12, 4]
    ref.species        <- mtag[mtag[,2] == 18, 4]
    particle.no        <- mtag[mtag[,2] == 13, 4]
    bins               <- mtag[mtag[,2] == 16, 4]
    
    seek(to.read, where = mtag[mtag[,2] == 5,1] + 8, origin = 'start')
    projectile         <- rawToChar(readBin(	to.read, raw(), n = mtag[mtag[,2] == 5,3], endian = endian))	
    seek(to.read, where = mtag[mtag[,2] == 4,1] + 8, origin = 'start')
    target.material    <- rawToChar(readBin(	to.read, raw(), n = mtag[mtag[,2] == 4,3], endian = endian))	
    seek(to.read, where = mtag[mtag[,2] == 6,1] + 8, origin = 'start')
    beam.energy.MeV.u  <- readBin( to.read, double(), size = 8, n = floor(mtag[mtag[,2] == 6,3]/8), endian = endian)
    seek(to.read, where = mtag[mtag[,2] == 7,1] + 8, origin = 'start')
    peak.position.g.cm2<- readBin( to.read, double(), size = 8, n = floor(mtag[mtag[,2] == 7,3]/8), endian = endian)
    
    if(header.only == FALSE){
      mm                 <- matrix(data = 0, ncol = 10, nrow = sum(bins))
      mm[,1]             <- rep.int(rep.int(1:n.depth.steps, n.particle.species), bins)
      mm[,2]             <- rep.int(rep.int(depth.g.cm2, n.particle.species), bins)
      mm[,3]             <- rep.int(sequence(n.particle.species), bins)
      mm[,4]             <- rep.int(particle.no, bins)
      
      fluence.tags       <- mtag[mtag[,2] == 19,]
      idx                <- 1
      for(i in 1:nrow(fluence.tags)){
        # i <- 1
        seek(to.read, where = fluence.tags[i,1] + 8, origin = 'start')
        size                   <- floor(fluence.tags[i,3]/8)
        mm[idx:(idx+size-1),9] <- readBin( to.read, double(), size = 8, n = size, endian = endian)
        idx                    <- idx + size
      }
      
      E.grid.tags        <- cbind( mtag[mtag[,2] %in% c(17,18),], 
                                   rep.int(1:n.depth.steps, n.particle.species),   # depth step
                                   sequence(n.particle.species),                   # species
                                   cumsum(bins)-bins[1]+1,                         # index in mm
                                   bins)                                           # size
      idx                <- 1
      for(i in 1:nrow(E.grid.tags)){
        # i <- 2
        if(E.grid.tags[i,2] == 17){
          seek(to.read, where = E.grid.tags[i,1] + 8, origin = 'start')
          size                   <- floor(E.grid.tags[i,3]/8)
          E.bins.MeV.u           <- readBin( to.read, double(), size = 8, n = size, endian = endian)
          E.low.MeV.u            <- E.bins.MeV.u[-length(E.bins.MeV.u)]
          E.high.MeV.u           <- E.bins.MeV.u[-1]
          mm[idx:(idx+size-2),5] <- E.low.MeV.u
          if(mean == "geometric"){  
            mm[idx:(idx+size-2),6]  <- sqrt(E.low.MeV.u * E.high.MeV.u)
          }else{
            mm[idx:(idx+size-2),6]  <- (E.low.MeV.u + E.high.MeV.u)/2
          }            
          mm[idx:(idx+size-2),7] <- E.high.MeV.u
          mm[idx:(idx+size-2),8] <- E.high.MeV.u - E.low.MeV.u
          idx                    <- idx + size - 1
        }else{
          ii                      <- E.grid.tags[,5] == E.grid.tags[i,5]
          ref.idx                 <- which((E.grid.tags[i,4]+1) == E.grid.tags[ii,6])
          size                    <- E.grid.tags[ii,8][ref.idx]
          from                    <- E.grid.tags[ii,7][ref.idx]
          mm[idx:(idx+size-1),5]  <- mm[from:(from+size-1),5]
          mm[idx:(idx+size-1),6]  <- mm[from:(from+size-1),6]
          mm[idx:(idx+size-1),7]  <- mm[from:(from+size-1),7]
          mm[idx:(idx+size-1),8]  <- mm[from:(from+size-1),8]
          idx                     <- idx + size
        }
      }
      mm[,10]      <- mm[,9] * mm[,8]         # convert fluence / binwidth -> fluence
      mm           <- mm[,-3]
      df           <- as.data.frame(mm)
      names(df)    <- c("depth.step", "depth.g.cm2", "particle.no", "E.low.MeV.u", "E.MeV.u", "E.high.MeV.u", "dE.MeV.u", "dN.dE.per.MeV.u.per.primary", "N.per.primary")
      cat(paste("Read ", n.depth.steps, " depth steps for projectile ", projectile, " on ", target.material, " with ", beam.energy.MeV.u, " MeV/u and peak at ", peak.position.g.cm2, " g/cm2.\n", sep = ""))
    }else{
      df           <- NULL
    }
    
    close(to.read)
    
  }else{
    # C version
    E.MeV.u             <- numeric(1)
    peak.position.g.cm2 <- numeric(1)
    particle.no         <- integer(1)
    material.no         <- integer(1)
    normalization       <- numeric(1)
    n.depth.steps       <- integer(1)
    returnValue         <- integer(1)
    res                 <- .C( "AT_SPC_read_header_from_filename_fast_R",
                               file.name           = as.character(file.name),
                               E.MeV.u             = as.single(E.MeV.u),
                               peak.position.g.cm2 = as.single(peak.position.g.cm2),
                               particle.no         = as.integer(particle.no),
                               material.no         = as.integer(material.no),
                               normalization       = as.single(normalization),
                               n.depth.steps       = as.integer(n.depth.steps),
                               returnValue         = as.integer(returnValue),
                               PACKAGE             = "libamtrack")
    beam.energy.MeV.u   <- res$E.MeV.u
    peak.position.g.cm2 <- res$peak.position.g.cm2
    n.depth.steps       <- res$n.depth.steps
    
    # TODO: projectile/material is hardcoded, replace by more flexible code
    projectile          <- "12C"
    target.material     <- "H2O"
    
    if(header.only == TRUE){
      df                  <- 0
    }else{
      spc.size            <- numeric(1)
      res                 <- .C( "AT_SPC_get_number_of_bins_from_filename_fast_R",
                                 file.name          = as.character(file.name),
                                 spc.size           = as.integer(spc.size),
                                 PACKAGE            = "libamtrack")
      
      n                    <- res$spc.size
      depth.step           <- integer(n)
      depth.g.cm2          <- numeric(n)
      E.MeV.u              <- numeric(n)
      DE.MeV.u             <- numeric(n)
      particle.no          <- integer(n)
      fluence.cm2          <- numeric(n)
      n.bins.read          <- integer(1)
      
      res                  <- .C( "AT_SPC_read_data_from_filename_fast_R",
                                  file.name          = as.character(file.name),
                                  n                  = as.integer(n),
                                  depth.step         = as.integer(depth.step),
                                  depth.g.cm2        = as.single(depth.g.cm2),
                                  E.MeV.u            = as.single(E.MeV.u),
                                  DE.MeV.u           = as.single(DE.MeV.u),
                                  particle.no        = as.integer(particle.no),
                                  fluence.cm2        = as.single(fluence.cm2),
                                  n.bins.read        = as.integer(n.bins.read),
                                  PACKAGE            = "libamtrack")
      
      df   <- data.frame(  depth.step             = res$depth.step,
                           depth.g.cm2            = res$depth.g.cm2,
                           E.MeV.u                = res$E.MeV.u,
                           DE.MeV.u               = res$DE.MeV.u,
                           particle.no            = res$particle.no,
                           fluence.cm2            = res$fluence.cm2)
      
    }
  }
  return(list( spc                 = df,
               n.depth.steps       = n.depth.steps,
               projectile          = projectile,
               target.material     = target.material,
               energy.MeV.u        = beam.energy.MeV.u,
               peak.position.g.cm2 = peak.position.g.cm2))
}

#############################
# tapply version for spectra
SPC.tapply <- function( x, 
                        INDEX, 
                        FUN, 
                        mixed.field.arguments = list(E.MeV.u     = "E.mid.MeV.u", 
                                                     fluence.cm2 = "N.per.primary", 
                                                     particle.no = "particle.no"),
                        additional.arguments = NULL, 
                        names.results        = NULL)
{    
  # Get index columns and levels
  index.columns    <- which(is.element(names(x@spectra), INDEX))
  if(length(INDEX) != length(index.columns)){
    cat("At least one index variable not found in spc data.\n")
    return(NULL)
  }
  
  index.variable   <- NULL
  for (i in 1:length(index.columns)){
    # DEBUG: i <- 1
    index.variable    <- paste(index.variable, x@spectra[,index.columns[i]])
  }
  levels           <- unique(index.variable)
  
  ###########################
  # Get argument list for FUN
  # Match with mixed field args
  args.FUN         <- names(formals(FUN))
  args.list        <- "("
  for (i in 1:length(args.FUN)) {
    mixed.field.arguments.idx <- match(args.FUN[i], names(mixed.field.arguments))
    if (!is.na(mixed.field.arguments.idx)) {
      args.list <- paste(	args.list, 
                          args.FUN[i], " = x@spectra$", 
                          mixed.field.arguments[[mixed.field.arguments.idx]], "[ii],", 
                          sep = "")
    }
  }
  
  if(!is.null(additional.arguments)){
    for(j in 1:length(additional.arguments)){
      if(additional.arguments[[j]][3] == TRUE){
        args.list    <- paste( args.list, 
                               additional.arguments[[j]][1], 
                               " = x@spectra$",
                               additional.arguments[[j]][2],
                               "[ii],",
                               sep = "")
      }else{
        args.list    <- paste( args.list, 
                               additional.arguments[[j]][1], 
                               " = ",
                               additional.arguments[[j]][2],
                               ",",
                               sep = "")
      }
    }
  }
  args.list        <- paste(substring(args.list, 1, nchar(args.list) - 1), ")")
  
  df.return        <- NULL
  for(cur.level in levels){
    # DEBUG: cur.level <- levels[2]
    ii            <- index.variable == cur.level
    
    res           <- eval( parse( text = paste( "FUN",
                                                args.list,
                                                sep = "")))
    df.cur.level  <-  cbind.data.frame( unique(data.frame( x@spectra[ii,index.columns])), 
                                        res)
    if(is.null(df.return)){
      df.return    <- df.cur.level
    }else{
      df.return    <- rbind.data.frame( df.return,
                                        df.cur.level)
    }
  }
  row.names(df.return) <- 1:nrow(df.return)
  names(df.return)[1:length(index.columns)]  <- INDEX
  if(!is.null(names.results)){
    names(df.return)     <- c(names(df.return)[1:length(index.columns)], names.results)
  }
  return(df.return)
}


redistribute.spc <- function(spc, E.min.MeV.u = 0, E.max.MeV.u = 600, dE.MeV.u = 1.0){
  if(spc@redistributed == TRUE){
    warning("spc is already redistributed.")
    return(spc)
  }
  
  particle.nos <- sort(unique(spc@spectra[,"particle.no"]))
  if(any(!(particle.nos%in%c(1002,2004,3006,4008,5010,6012)))){
    stop("SPC contains more than six canonical particles.")
  }
  
  depth.steps  <- sort(unique(spc@spectra[,"depth.step"]))
  
  E.MeV.u      <- seq(E.min.MeV.u + dE.MeV.u/2, 
                      E.max.MeV.u - dE.MeV.u/2,
                      dE.MeV.u)
  
  cv                        <- paste(spc@spectra[,"depth.step"], spc@spectra[,"particle.no"])
  n.E                       <- length(E.MeV.u)
  
  new.spectra               <- spectra(n.E * length(particle.nos) * length(depth.steps))
  
  new.spectra[,"depth.step"]<- unlist( tapply( spc@spectra[,"depth.step"],
                                               cv,
                                               function(x, n){rep(unique(x), n)},
                                               n = n.E))
  new.spectra[,"depth.g.cm2"] <- unlist( tapply( spc@spectra[,"depth.g.cm2"],
                                                                         cv,
                                                                         function(x, n){rep(unique(x), n)},
                                                                         n = n.E))
  new.spectra[,"particle.no"] <-  unlist( tapply( spc@spectra[,"particle.no"],
                                                                         cv,
                                                                         function(x, n){rep(unique(x), n)},
                                                                         n = n.E))
  new.spectra[,"E.MeV.u"]     <-  unlist( tapply( spc@spectra[,"E.MeV.u"],
                                                                         cv,
                                                                         function(x, E){E},
                                                                         E = E.MeV.u))
  new.spectra[,"dE.MeV.u"]    <- dE.MeV.u
  new.spectra[,"N.per.primary"] <-  unlist( by(   spc@spectra,
                                                                         cv,
                                                                         function(x, E, dE){
                                                                           yy <- approx(x    = x[,"E.MeV.u"],
                                                                                        y    = x[,"N.per.primary"] / x[,"dE.MeV.u"],
                                                                                        xout = E)$y * dE
                                                                           yy[is.na(yy)] <- 0.0
                                                                           yy},
                                                                         E = E.MeV.u,
                                                                         dE = dE.MeV.u))
  cat("spc redistributed with E.min", E.min.MeV.u, "MeV/u, E.max", E.max.MeV.u, "MeV/u and dE", dE.MeV.u, "\n")
  return( new( "dataSPC",
               projectile          = spc@projectile,
               beam.energy.MeV.u   = spc@beam.energy.MeV.u,
               target.material     = spc@target.material,
               peak.position.g.cm2 = spc@peak.position.g.cm2,
               redistributed       = TRUE, 
               spectra             = new.spectra))
}

