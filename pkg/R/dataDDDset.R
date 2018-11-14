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
dataDDDset <- function(type = c("ddd", "csv")[1], ddd.path = "."){
    
  if(type == "ddd"){
    pattern <- "*\\.ddd"
  }
  if(type == "csv"){
    pattern <- "*\\.csv"
  }

  files <- list.files(path       = ddd.path, 
                      pattern    = pattern,
                      full.names = TRUE)
  if(length(files) == 0){
    stop("No ddd files found - wrong directory?")
  }

  DDDs                   <- lapply(files, function(x){
                                            # x <- files[1]
                                            cat("Reading", basename(x), "...\n")
                                            dataDDD(x)})
  
  energies               <- sapply(DDDs, function(x){x@beam.energy.MeV.u})
  DDDs                   <- DDDs[order(energies)]
  
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
  
  if(DDD.set@beam.energies.MeV.u[closest.idx] != beam.energy.MeV.u){
    cat("Closest energy available is", DDD.set@beam.energies.MeV.u[closest.idx], "MeV/u, returning this DDD.\n")}
  
  return(DDD.set@DDDs[[closest.idx]])
}


get.dose.Gy.from.set <- function(DDD.set, depths.g.cm2, weights){
  
  n.ddd <- length(DDD.set@beam.energies.MeV.u)
  
  if (missing(weights)){
    weights <- rep(1.0, n.ddd)
  }
  
  doses     <- matrix( unlist( mclapply( 1:n.ddd, 
                                       function(i, x, y){ get.dose.Gy( get.ddd( x, 
                                                                                x@beam.energies.MeV.u[i]), 
                                                                        y)},
                                       x = DDD.set,
                                       y = depths.g.cm2)),
                       nrow  = n.ddd,
                       byrow = TRUE)
  return(apply(weights * doses, 2, sum))
}

setMethod(f          = "[", 
          signature  = "dataDDDset",
          definition = function(x, i){
            
            # Check if binary image is given for masking
            if(!missing(i)){
              return(new("dataDDDset",
                               projectiles          = x@projectiles[i],
                               beam.energies.MeV.u  = x@beam.energies.MeV.u[i],
                               target.materials     = x@target.materials[i],
                               peak.positions.g.cm2 = x@peak.positions.g.cm2[i],
                               densities.g.cm3      = x@densities.g.cm3[i],
                               DDDs                 = x@DDDs[i]))
            }else{
                return(x)
              }
            }
)

          