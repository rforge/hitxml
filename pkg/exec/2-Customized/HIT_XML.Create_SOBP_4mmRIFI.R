#' Script to generate a spread-out Bragg peak
#'
#' Optimized either wrt physical or biological dose. 
#' 
#' Original version by Grischa Klimpki, 2014
#' adapted to HIT_XML package and extended to 
#' biological optization by Steffen Greilich, 2015-07
#' 
#' Revised 2017-09-25

rm(list = ls())
library(HITXML)
library(lattice)
library(parallel)


#' START OF USER INPUT

# path to spc and ddd data
#ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/p/RF0MM/"
#spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/p/RF0MM/"
ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/06 - RIFI 4mm/ddd_12C-RF4mm/"
spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM/"
rbe.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/RBE"

# minimal and maximal depth in cm
min.depth.g.cm2         <- 10
max.depth.g.cm2         <- 16
expid                   <- "SOBP_4mmRIFI"

offset.g.cm2            <- 0.289

step.size.g.cm2         <- 0.025
IES.step                <- 5
plateau.dose.Gy         <- 1

output.LET              <- FALSE
LET.step.size.g.cm2     <- 0.125

biol.optimization       <- FALSE
rbe.file                <- "chordom02.rbe"
n.biol.opt.steps        <- 5
bio.step.size.g.cm2     <- 0.25

write.SOBP              <- TRUE

plot.range              <- 2.0

#' END OF USER INPUT

# Add offset to SOBP range
min.depth.g.cm2         <- min.depth.g.cm2 + offset.g.cm2
max.depth.g.cm2         <- max.depth.g.cm2 + offset.g.cm2

# read in DDD files  
ddds              <- dataDDDset(ddd.path = ddd.path)

# get IESs necessary
jj                <- ddds@peak.positions.g.cm2 > min.depth.g.cm2 & ddds@peak.positions.g.cm2 <= max.depth.g.cm2 &
                             rep(c(TRUE, rep(FALSE, IES.step - 1)), length.out = length(ddds@projectiles))
#jj[head(which(jj), 1)-IES.step] <- TRUE
#jj[tail(which(jj), 1)+IES.step] <- TRUE
no.IES            <- sum(jj)
ddds.sub          <- ddds[which(jj)]


#' A. DO PHYSICAL OPTIMIZATION

# Vector of depths covering the SOBP
depths.g.cm2       <- seq( from       = min.depth.g.cm2, 
                           to         = max.depth.g.cm2, 
                           by         = step.size.g.cm2)

# Get plateau per primary
plateau.dose.per.primary.Gy <- mean(get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                                         depths.g.cm2 = depths.g.cm2))

# Objective function: sum of squares (or higher) of deviation to dose set
# p <- rep( 1, no.IES) + seq(0, 1, length.out = no.IES)
# DDD.set <- ddds.sub
# dose.set.Gy <- plateau.dose.per.primary.Gy

dose.dev <- function(p, depths.g.cm2, DDD.set, dose.set.Gy){
  doses <- get.dose.Gy.from.set(DDD.set      = DDD.set, 
                                depths.g.cm2 = depths.g.cm2, 
                                weights      = p)
  idx <- 1:length(depths.g.cm2)
  
#  weights.of.weights <- (idx - mean(idx))^2
#  cost <- log10(sum( weights.of.weights * (doses - dose.set.Gy)^2))
  cost <- log10(sum((doses - dose.set.Gy)^2))
  return(cost)
}

# minimize objective function to find weights
rel.weights       <- optim( fn             = dose.dev, 
                            par            = rep( 1, no.IES) + seq(0, 1, length.out = no.IES),
                            depths.g.cm2   = depths.g.cm2,
                            DDD.set        = ddds.sub,
                            dose.set.Gy    = plateau.dose.per.primary.Gy,
                            method         = "L-BFGS-B",
                            lower          = rep(0.1, no.IES),
                            upper          = rep(20, no.IES),
                            control        = list(trace = TRUE, 
                                                  maxit = 300,
                                                  factr = 1e9))$par

# Scale number of primaries to get actual weights
fluence.factor <- plateau.dose.Gy / plateau.dose.per.primary.Gy
total.weights  <- rel.weights * fluence.factor

##########################
# plot SOBP (single field) 
plot.depths.g.cm2       <- seq( from       = 0.0, 
                                to         = max(ddds.sub@peak.positions.g.cm2) * plot.range, 
                                by         = step.size.g.cm2)

print(plot.SOBP(plot.ddds = ddds.sub, 
               plot.depths.g.cm2 = plot.depths.g.cm2, 
               plot.weights = total.weights, 
               "(phys.opt.)",
               start.depth.cm = min.depth.g.cm2,
               end.depth.cm = max.depth.g.cm2))


#################################
# B. Biological optimzation steps
if(biol.optimization | output.LET){
  #################
  # 1. Read in data
  
  # Scan available spc files to get energy sample points
  spcs             <- dataSPCset(spc.path = spc.path)
  
  # Returns a list of dataSPC objects with same energies
  # than selected IESs. dataSPCs will be interpolated
  # between the adjacent energies of the available
  # energy sample point in the spc set
  ddds.spc          <- mclapply(ddds.sub@beam.energies.MeV.u, 
                                function(x, s){
                                  get.spc(s, x)},
                                spcs)
  # save(ddds.spc, file = "ddds.spc.rda")
  # load("ddds.spc.rda")
  
  rbe.data <- dataRBE(rbe.file, rbe.path)
  
  # no.IES <- 1
  # depths.g.cm2 <- 0
  

  if(output.LET){
    LET.depths.g.cm2 <- seq(0, max.depth.g.cm2*2, by = LET.step.size.g.cm2)
    
    # Returns a list (length number of IESs in SOBP)
    # of lists (each length number of depths covering entire depth dose curve)
    # of dataSpectrum objects - representing the spectra
    # from the individual IESs (first index) contributing at specific depth (second index)
    # TODO: Move this format to its own class and adapt lapply's below
    spectra.at.depth  <- mclapply(1:no.IES,
                                  function(i, s, d){
                                    cat("Getting spectra from IES", i, "\n")
                                    dataSpectrum(s[[i]], d)},
                                  s = ddds.spc,
                                  d = LET.depths.g.cm2)
    # save(spectra.at.depth, file = "sad.rda")
    # load("sad.rda")
    
    w.spectra.at.depth  <- mclapply(1:no.IES,
                                    function(i, s, w){
                                      cat("Weighting spectra from IES", i, "\n")
                                      lapply(1:length(s[[i]]),
                                           function(j, ss, ww){
                                             ss[[j]] * ww},
                                           ss = s[[i]],
                                           ww = w[i])},
                                  s = spectra.at.depth,
                                  w = total.weights)
    
    eff.spectra.at.depth <- mclapply(1:length(LET.depths.g.cm2),
                                   function(i, s, n){
                                     cat("Adding spectra for depth", i, "\n")
                                     ss <- s[[1]][[i]]
                                     if(n > 1){
                                       for(j in 2:n){
                                         ss <- ss + s[[j]][[i]]
                                       }
                                     }
                                     ss},
                                   s = w.spectra.at.depth,
                                   n = no.IES)
    
    total                       <- sapply(eff.spectra.at.depth,
                                          spectrum.total.n.particles)
    primaries                   <- sapply(eff.spectra.at.depth,
                                          spectrum.total.n.particles,
                                          particle.no = 6012)
    Z1                          <- sapply(eff.spectra.at.depth,
                                          spectrum.total.n.particles,
                                          particle.no = 1002)
    Z2                          <- sapply(eff.spectra.at.depth,
                                          spectrum.total.n.particles,
                                          particle.no = 2004)
    fLET.total                  <- sapply(eff.spectra.at.depth,
                                          spectrum.fLET)/10
    fLET.primaries              <- sapply(eff.spectra.at.depth,
                                          spectrum.fLET,
                                          particle.no = 6012)/10
    dLET.total                  <- sapply(eff.spectra.at.depth,
                                          spectrum.dLET)/10

    fLET.primaries[LET.depths.g.cm2 > max.depth.g.cm2] <- NA

    D.phys.Gy                  <-  get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                                        depths.g.cm2 = LET.depths.g.cm2, 
                                                        weights      = total.weights)

    rbe                   <- HX.RBE.LEM(rbe.data, 
                                        eff.spectra.at.depth, 
                                        D.phys.Gy)$RBE

    rbe[LET.depths.g.cm2 > max.depth.g.cm2 *1.5] <- NA
    
    xyplot(total + primaries + Z1 + Z2 ~ LET.depths.g.cm2,
           grid = TRUE,
           type = "l",
           auto.key = list(space = "right"),
           xlab = list("depth / cm", cex=1.5),
           ylab = list("fluence / (1/cm2)", cex=1.5),
           scale = list(cex = 1.25))

    xyplot(fLET.primaries + fLET.total + dLET.total ~ LET.depths.g.cm2,
           grid = TRUE,
           type = "l",
           auto.key = list(space = "right"),
           xlab = list("depth / cm", cex=1.5),
           ylab = list("fluence / (1/cm2)", cex=1.5),
           scale = list(cex = 1.25))
    
    xyplot(rbe ~ LET.depths.g.cm2,
           grid = TRUE,
           type = "l",
           auto.key = list(space = "right"),
           xlab = list("depth / cm", cex=1.5),
           ylab = list("RBE", cex=1.5),
           scale = list(cex = 1.25))
    
    D.bio.GyRBE <- D.phys.Gy * rbe
    
    xyplot(D.phys.Gy + D.bio.GyRBE ~ LET.depths.g.cm2,
           grid = TRUE,
           type = "l",
           auto.key = list(space = "right"),
           xlab = list("depth / cm", cex=1.5),
           ylab = list("dose / Gy, GyRBE", cex=1.5),
           scale = list(cex = 1.25))

    if(!biol.optimization){
      stop("No error! Just no biological optimization requested")
    }
  }

  # Returns a list (length number of IESs in SOBP)
  # of lists (each length number of depths covering SOBP)
  # of dataSpectrum objects - representing the spectra
  # from the individual IESs (first index) contributing at specific depth (second index)
  # TODO: Move this format to its own class and adapt lapply's below
  spectra.at.depth  <- mclapply(1:no.IES,
                              function(i, s, d){
                                cat("Getting spectra from IES", i, "\n")
                                dataSpectrum(s[[i]], d)},
                              s = ddds.spc,
                              d = depths.g.cm2)
  # save(spectra.at.depth, file = "sad.rda")
  # load("sad.rda")

  
  ###############################
  # 2. Iterate with RBE weighting
  for(biol.opt.step in 1:n.biol.opt.steps){
    # Applies the given weights to all spectra at depths covering SOBP (second index)
    # from respective IESs (first index) 
    w.spectra.at.depth  <- mclapply(1:no.IES,
                                  function(i, s, w){
                                    cat("Weighting spectra from IES", i, "\n")
                                    lapply(1:length(s[[i]]),
                                           function(j, ss, ww){
                                             ss[[j]] * ww},
                                           ss = s[[i]],
                                           ww = w[i])},
                                  s = spectra.at.depth,
                                  w = total.weights)
    
    # Combines (adds) spectra from all IESs (first index) for the depths 
    # covering the SOBP (second index). Spectra can be have been weighted before
    eff.spectra.at.depth <- mclapply(1:length(depths.g.cm2),
                 function(i, s, n){
                   cat("Adding spectra for depth", i, "\n")
                   ss <- s[[1]][[i]]
                   if(n > 1){
                     for(j in 2:n){
                      ss <- ss + s[[j]][[i]]
                     }
                   }
                   ss},
                 s = w.spectra.at.depth,
                 n = no.IES)
    
    # Why does this get different results??
    #dose.Gy(pp, 
    #        stopping.power.source = "ICRU", 
    #        target.material = "Water, Liquid")
    D.phys.Gy             <- get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                   depths.g.cm2 = depths.g.cm2, 
                                                   weights      = total.weights)
    rbe                   <- HX.RBE.LEM(rbe.data, 
                                        eff.spectra.at.depth, 
                                        D.phys.Gy)$RBE
    
    # Get current physical dose in plateau per Gy
    plateau.dose.per.primary.Gy <- mean(get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                                             depths.g.cm2 = depths.g.cm2))
    # Scale down using RBE at all depth of SOBP
    rbe.scaled.plateau.doses.per.primary.Gy <- plateau.dose.per.primary.Gy / rbe

    # minimize objective function to find weights
    rel.weights       <- optim( fn             = dose.dev, 
                                par            = rep( 1, no.IES),
                                depths.g.cm2   = depths.g.cm2,
                                DDD.set        = ddds.sub,
                                dose.set.Gy    = rbe.scaled.plateau.doses.per.primary.Gy,
                                method         = "L-BFGS-B",
                                lower          = rep(0.0, no.IES),
                                control        = list(trace = TRUE, 
                                                      maxit = 200,
                                                      reltol = 1e-3))$par
    # Scale number of primaries to get actual weights
    fluence.factor <- plateau.dose.Gy / mean(plateau.dose.per.primary.Gy)
    total.weights  <- rel.weights * fluence.factor
    
    plot(plot.SOBP(ddds.sub, 
                   plot.depths.g.cm2, 
                   total.weights, 
                   paste0("(biol.opt., step ", biol.opt.step, ")")))
  }

  ###########
  # GET RBE for entire field ( with reduced resolution)
  bio.depths.g.cm2         <- seq( from       = 0.0, 
                                  to         = max(ddds.sub@peak.positions.g.cm2) * plot.range, 
                                  by         = bio.step.size.g.cm2)
  
  spectra.at.depth  <- mclapply(1:no.IES,
                              function(i, s, d){
                                cat("Getting spectra from IES", i, "\n")
                                dataSpectrum(s[[i]], d)},
                              s = ddds.spc,
                              d = bio.depths.g.cm2)

  w.spectra.at.depth  <- mclapply(1:no.IES,
                                function(i, s, w){
                                  cat("Weighting spectra from IES", i, "\n")
                                  lapply(1:length(s[[i]]),
                                         function(j, ss, ww){
                                           ss[[j]] * ww},
                                         ss = s[[i]],
                                         ww = w[i])},
                                s = spectra.at.depth,
                                w = total.weights)
  
  eff.spectra.at.depth <- mclapply(1:length(bio.depths.g.cm2),
                                 function(i, s, n){
                                   cat("Adding spectra for depth", i, "\n")
                                   ss <- s[[1]][[i]]
                                   if(n > 1){
                                     for(j in 2:n){
                                       ss <- ss + s[[j]][[i]]
                                     }
                                   }
                                   ss},
                                 s = w.spectra.at.depth,
                                 n = no.IES)
  
  D.phys.Gy             <- get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                 depths.g.cm2 = bio.depths.g.cm2, 
                                                 weights      = total.weights)
  
  rbe                   <- HX.RBE.LEM(rbe.data, eff.spectra.at.depth, D.phys.Gy)$RBE

  df.plot <- data.frame( depth.g.cm2          = plot.depths.g.cm2,
                         D.phys.Gy            = get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                                      depths.g.cm2 = plot.depths.g.cm2, 
                                                                      weights      = total.weights),
                         RBE                  = approx( x    = bio.depths.g.cm2,
                                                        y    = rbe,
                                                        xout = plot.depths.g.cm2,
                                                        rule = 2)$y)
  df.plot$D.biol.Gy <- df.plot$D.phys.Gy * df.plot$RBE
  
  # myyscale.component <- function(...) 
  # { 
  #   ans 				              <- yscale.components.default(...) 
  #   ans$right 	              <- ans$left 
  #   foo 				              <- ans$right$labels$at 
  #   ans$right$labels$labels 	<- as.character(foo * 5) 
  #   ans 
  # } 
  # 
  # xyplot(D.biol.Gy + D.phys.Gy + RBE/5 ~ depth.g.cm2,
  #        df.plot,
  #        type = "l",
  #        xlab = list("depth / (g/cm2)", cex=1.5),
  #        ylab = list("Dose / (Gy, GyRBE)", cex=1.5),
  #        par.settings=	list( layout.widths=list(right.padding=10)), 
  #        yscale.component	= myyscale.component,
  #        legend = list(right = list(fun = textGrob, 
  #                                   args = list(x = 3, 
  #                                               y = 0.5, 
  #                                               rot = 90,
  #                                               label = "RBE",
  #                                               gp=gpar(cex = 1.5, col = "darkgreen"),
  #                                               just = "center", 
  #                                               default.units = "native", 
  #                                               vp = viewport(xscale = c(0, 1), 
  #                                                             yscale = c(0, 1))))),
  #        scales = list( x = list(cex = 1.25),
  #                       y = list(cex = 1.25, 
  #                                relation = "free", 
  #                                rot = 0)),
  #        grid = TRUE,
  #        main = list(paste0("SOBP (single field, ",
  #                           unique(ddds.sub@projectiles),
  #                           ") consisting of ", 
  #                           no.IES, 
  #                           " IESs "), 
  #                    cex=1.5),
  #        sub  = "NB: HIT isocenter is at 0.289 g/cm2 depth!")
}

###############
# save SOBP.dat
if(write.SOBP){
  df                 <- data.frame(E.GeV   = ddds.sub@beam.energies.MeV.u / 1000,
                                   x.cm    = rep(0.0, no.IES),
                                   y.cm    = rep(0.0, no.IES),
                                   FWHM.cm = rep(0.1, no.IES),
                                   fluence = as.integer(total.weights))
  file.name <- paste0("SOBP_", expid, "_simple.dat")
  write.table(df,
              file = file.name,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              eol       = "\r\n" )
  cat("Resulting weights written to ", file.name, "\n")
}

alarm()