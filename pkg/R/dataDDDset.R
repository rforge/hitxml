################################
# dataDDDset CLASS
################################
setClass( Class            = "dataDDDset",
          slots            = c( projectiles              = "character",
                                beam.energies.MeV.u      = "numeric",
                                target.materials         = "character",
                                peak.positions.g.cm2     = "numeric",
                                densities.g.cm3          = "numeric",
                                alpha.X.Gy               = "numeric",
                                beta.X.Gy2               = "numeric",
                                RBE.model                = "character",
                                DDDs                     = "list"),
          prototype        = list( projectiles            = character(),
                                   beam.energies.MeV.u    = numeric(),
                                   target.materials       = character(),
                                   peak.positions.g.cm2   = numeric(),
                                   densities.g.cm3        = numeric(),
                                   alpha.X.Gy             = numeric(),
                                   beta.X.Gy2             = numeric(),
                                   RBE.model              = character(),
                                   DDDs                   = list()) )

################################
# Constructor
dataDDDset <- function(ddd.path = "."){
    
  files <- list.files(path       = ddd.path, 
                      pattern    = "*\\.ddd",
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
  
  alpha.X.Gy             <- unique(sapply(DDDs, function(x){x@alpha.X.Gy}))
  beta.X.Gy2             <- unique(sapply(DDDs, function(x){x@beta.X.Gy2}))
  RBE.model              <- unique(sapply(DDDs, function(x){x@RBE.model}))
  if(length(alpha.X.Gy) > 1){
    stop("DDD data differ in alphaX value!")
  }
  if(length(beta.X.Gy2) > 1){
    stop("DDD data differ in betaX value!")
  }
  if(length(RBE.model) > 1){
    stop("DDD data differ in RBE model!")
  }

  new("dataDDDset",
      projectiles          = sapply(DDDs, function(x){x@projectile}),
      beam.energies.MeV.u  = sapply(DDDs, function(x){x@beam.energy.MeV.u}),
      target.materials     = sapply(DDDs, function(x){x@target.material}),
      peak.positions.g.cm2 = sapply(DDDs, function(x){x@peak.position.g.cm2}),
      densities.g.cm3      = sapply(DDDs, function(x){x@density.g.cm3}),
      alpha.X.Gy           = alpha.X.Gy,
      beta.X.Gy2           = beta.X.Gy2,
      RBE.model            = RBE.model,
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


get.RBE.from.set <- function(DDD.set, depths.g.cm2, weights){
  
  n.ddd <- length(DDD.set@beam.energies.MeV.u)
  
  if (missing(weights)){
    weights <- rep(1.0, n.ddd)
  }
  
  # Get doses, alphas and betas for all weights (rows) and depths (cols)
  doses     <- matrix( unlist( mclapply( 1:n.ddd, 
                                         function(i, x, y){ get.dose.Gy( get.ddd( x, 
                                                                                  x@beam.energies.MeV.u[i]), 
                                                                         y)},
                                         x = DDD.set,
                                         y = depths.g.cm2)),
                       nrow  = n.ddd,
                       byrow = TRUE)

  alphas   <- matrix( unlist( mclapply( 1:n.ddd, 
                                         function(i, x, y){ 
                                           # i <- 1
                                           get.alpha.ion.Gy( get.ddd( x, 
                                                                      x@beam.energies.MeV.u[i]), 
                                                                         y)},
                                         x = DDD.set,
                                         y = depths.g.cm2)),
                       nrow  = n.ddd,
                       byrow = TRUE)
  
  
  betas   <- matrix( unlist( mclapply( 1:n.ddd, 
                                        function(i, x, y){ 
                                          # i <- 1
                                          get.beta.ion.Gy2( get.ddd( x, 
                                                                     x@beam.energies.MeV.u[i]), 
                                                            y)},
                                        x = DDD.set,
                                        y = depths.g.cm2)),
                      nrow  = n.ddd,
                      byrow = TRUE)
  
  
  # Get total doses for all depths and relative dose contributions for all weights (rows) and depths (cols)
  weighted.doses     <- weights * doses
  total.doses        <- apply(weights * doses, 2, sum)
  rel.weighted.doses <- t(t(weighted.doses) / total.doses)
  
  # Sanity check: sum over weights (i.e. has to be 1.0)
  # apply(rel.weighted.doses, 2, sum) - 1.0
  
  # Use Keller&Rossi additivity rule to get mixed field alpha and beta for all depths
  weighted.alphas    <- rel.weighted.doses * alphas
  weighted.betas     <- (rel.weighted.doses * sqrt(betas))^2
  total.ion.alphas   <- apply(weighted.alphas, 2, sum)
  total.ion.betas    <- apply(weighted.betas, 2, sum)
  
  # Get log survival fractions for ions
  logSFs.ion <- total.doses * total.ion.alphas + total.doses^2 * total.ion.betas
  
  # Get full survival curve for X-rays
  Xray.dose.Gy       <- seq(0, 30, length.out = 3000)
  logSF.Xray         <- Xray.dose.Gy * DDD.set@alpha.X.Gy + Xray.dose.Gy^2 * DDD.set@beta.X.Gy2

  # Find iso-effective X-ray doses
  Xray.iso.doses.Gy  <- approx(x = logSF.Xray, y = Xray.dose.Gy, xout = logSFs.ion)$y
    
  # Get RBE and return
  return(Xray.iso.doses.Gy / total.doses)
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
                               alpha.X.Gy           = x@alpha.X.Gy,
                               beta.X.Gy2           = x@beta.X.Gy2,
                               RBE.model            = x@RBE.model,
                               DDDs                 = x@DDDs[i]))
            }else{
                return(x)
              }
            }
)

          