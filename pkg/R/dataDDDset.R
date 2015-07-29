################################
# dataDDDset CLASS
################################
setClass( Class            = "dataDDDset",
          slots            = c( projectiles              = "character",
                                beam.energies.MeV.u      = "numeric",
                                target.materials         = "character",
                                peak.positions.g.cm2     = "numeric",
                                densities.g.cm3          = "numeric",
                                DDDs                     = "list"),
          prototype        = list( projectiles            = character(),
                                   beam.energies.MeV.u    = numeric(),
                                   target.materials       = character(),
                                   peak.positions.g.cm2   = numeric(),
                                   densities.g.cm3        = numeric(),
                                   DDDs                   = list()) )

################################
# Constructor
dataDDDset <- function(pattern = "*.ddd", ddd.path = "."){

    files <- list.files(path       = ddd.path, 
                        pattern    = pattern,
                        full.names = TRUE)
    
    projectiles            <- character()
    beam.energies.MeV.u    <- numeric()
    target.materials       <- character()
    peak.positions.g.cm2   <- numeric()
    densities.g.cm3        <- numeric()
    DDDs                   <- lapply(files, function(x){ cat("Reading", basename(x), "...\n")
                                                         dataDDD(x)})
      
    new("dataDDDset",
        projectiles          = sapply(DDDs, function(x){x@projectile}),
        beam.energies.MeV.u  = sapply(DDDs, function(x){x@beam.energy.MeV.u}),
        target.materials     = sapply(DDDs, function(x){x@target.material}),
        peak.positions.g.cm2 = sapply(DDDs, function(x){x@peak.position.g.cm2}),
        densities.g.cm3      = sapply(DDDs, function(x){x@density.g.cm3}),
        DDDs                 = DDDs)
}

get.ddd <- function(DDD.set, beam.energy.MeV.u){

  if(beam.energy.MeV.u < min(DDD.set@beam.energies.MeV.u) | beam.energy.MeV.u > max(DDD.set@beam.energies.MeV.u)){
    stop("Requested energy outside of grid of available data.")
  }
  
  # find closest match
  distance    <- abs(DDD.set@beam.energies.MeV.u - beam.energy.MeV.u)
  closest.idx <- which.min(distance)
  
  cat("Closest energy available is", DDD.set@beam.energies.MeV.u[closest.idx], "MeV/u, returning this DDD.")
  return(DDD.set@DDDs[[closest.idx]])
}
