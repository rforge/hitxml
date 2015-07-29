################################
# dataSPCset CLASS
################################
setClass( Class            = "dataSPCset",
          slots            = c( projectiles              = "character",
                                beam.energies.MeV.u      = "numeric",
                                target.materials         = "character",
                                peak.positions.g.cm2     = "numeric",
                                endian                   = "character",
                                files                    = "character"),
          prototype        = list( projectiles            = character(),
                                   beam.energies.MeV.u    = numeric(),
                                   target.materials       = character(),
                                   peak.positions.g.cm2   = numeric(),
                                   endian                 = character(),
                                   files                  = character()) )

################################
# Constructor
dataSPCset <- function(pattern = "*.spc", spc.path = ".", endian = "little"){

    files <- list.files(path       = spc.path, 
                        pattern    = pattern,
                        full.names = TRUE)
    projectiles            <- character()
    beam.energies.MeV.u    <- numeric()
    target.materials       <- character()
    peak.positions.g.cm2   <- numeric()
    endian                 <- endian
    
    for (file in files){
      # file <- files[1]
      tmp    <- SPC.read(file.name   = file,
                         endian      = endian,
                         header.only = TRUE)
      projectiles               <- c(projectiles, tmp$projectile)
      beam.energies.MeV.u       <- c(beam.energies.MeV.u, tmp$energy.MeV.u)
      target.materials          <- c(target.materials, tmp$target.material)
      peak.positions.g.cm2      <- c(peak.positions.g.cm2, tmp$peak.position.g.cm2)
      cat("Read ", basename(file), ".\n")
    }

    new("dataSPCset",
        projectiles          = projectiles,
        beam.energies.MeV.u  = beam.energies.MeV.u,
        target.materials     = target.materials,
        peak.positions.g.cm2 = peak.positions.g.cm2,
        endian               = endian,
        files                = files)
}

get.spc <- function(SPC.set, beam.energy.MeV.u){

  if(beam.energy.MeV.u < min(SPC.set@beam.energies.MeV.u) | beam.energy.MeV.u > max(SPC.set@beam.energies.MeV.u)){
    stop("Requested energy outside of grid of available data.")
  }
  
  # find two closest matches
  distance    <- abs(SPC.set@beam.energies.MeV.u - beam.energy.MeV.u)
  closest.idx <- match(head(sort(distance), 2), distance)
  if(closest.idx[1] == closest.idx[2]){
    closest.idx[2] <- closest.idx[1] + 1
  }
  
  # read and interpolate SPCs
  spc.lower   <- dataSPC(file.name        = SPC.set@files[closest.idx[1]],
                         endian           = SPC.set@endian)
  spc.upper   <- dataSPC(file.name        = SPC.set@files[closest.idx[2]],
                         endian           = SPC.set@endian)
  
  spc         <- SPC.interpolate( spc.lower    = spc.lower,
                                  spc.upper    = spc.upper,
                                  energy.MeV.u = beam.energy.MeV.u)
  return(spc)
}

SPC.interpolate <- function(spc.lower, spc.upper, energy.MeV.u){
  if(length(unique(spc.lower@spectra$depth.step)) != length(unique(spc.upper@spectra$depth.step))){
    # TODO: rescale spc with more steps to steps from spc with
    # TODO: fewer steps using AT.SPC.spectrum.at.depth.g.cm2
    stop("Cannot interpolate, different number of depth steps")
  }
  if((spc.lower@projectile      != spc.upper@projectile)|
     (spc.lower@target.material != spc.upper@target.material)){
    stop("Cannot interpolate, target materials or projectiles do not match.")
  }

    # Copy structure
  spc                         <- spc.lower
  spc@beam.energy.MeV.u       <- energy.MeV.u
  
  # Get relative difference of requested energy
  frac                   <- (energy.MeV.u - spc.lower@beam.energy.MeV.u) / (spc.upper@beam.energy.MeV.u - spc.lower@beam.energy.MeV.u)
  
  # scale depth steps
  # TODO: should scales with E2 like range
  spc@spectra$depth.g.cm2    <- (1-frac) * spc.lower@spectra$depth.g.cm2 + frac * spc.upper@spectra$depth.g.cm2
  
  # interpolate fluences
  # TODO: check if linear interpolation really applies
  spc@spectra$N.per.primary    <- (1-frac) * spc.lower@spectra$N.per.primary + frac * spc.upper@spectra$N.per.primary
  
  return(spc)
}