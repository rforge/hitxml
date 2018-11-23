#' Script to generate a spread-out Bragg peak
#'
#' Optimized either wrt physical or biological dose. 
#' 
#' Original version by Grischa Klimpki, 2014
#' adapted to HIT_XML package and extended to 
#' biological optization by Steffen Greilich, 2015-07
#' 
#' Revised 2017-09-25
#' Revised 2018-06-23, write out of depth curves, spectra, SG
#' Revised 2018-11-13, use generic RBE data
rm(list = ls())
library(HITXML)
library(lattice)
library(parallel)
library(optimParallel)
library(data.table)
library(grid)

#===============#
# USER INPUT ####
#===============#
expid                   <- "sg86110"
# path to spc and ddd data
base.path <- file.path("D:/00 - Einstellungen/E0409-NB7/",
                       "Dropbox/Beruf/Workspace/TRS398_SPR_revision/03 - Data/TRS398_C12_basedata_generic/generic")

ddd.path <- file.path(base.path,
                      "DDD-HITXML/12C/RF3MM_3mmSteps_enSpread/ab2")

#spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM/"

# minimal and maximal depth in cm
min.depth.g.cm2         <- 6.0
max.depth.g.cm2         <- 16.0
extend.IES              <- c(1,1)   # Number of additional IESs proximal and distal to SOBP (can yield more smooth plateau) 

# WET offset to isocenter
offset.g.cm2            <- 0.0

# SOBP plateau
step.size.g.cm2         <- 0.02     # Distance between dose points which enter plateau flatness optimization 
IES.step                <- 1        # Use every ... IES available
plateau.dose.Gy         <- 3        # Dose at SOBP plateau (GyRBE)

output.LET              <- FALSE
LET.step.size.g.cm2     <- 0.125

biol.optimization       <- TRUE
rbe.method              <- c("SPCs and RBE file", "From alpha/beta with depth")[2]
rbe.path                <- file.path(base.path,
                                     "DDD_alpha_beta/12C/RF3MM_3mmSteps_enSpread")
rbe.file                <- "chordom02.rbe"
n.biol.opt.steps        <- 10
bio.step.size.g.cm2     <- 0.02

# Output spectra?
output.spectra          <- FALSE     
spectra.step.size.g.cm2 <- 0.125    # Distance between depth positions at which spectra are taken

# Misc
plot.range              <- 2.0      # Plot depth up to ... time Bragg peak position 
conversion.threshold    <- 1e-9     # Criterion of conversion for cost function (<1e-7)


#==========================#
# Physical optimization ####
#==========================#

# Add offset to SOBP range
min.depth.g.cm2         <- min.depth.g.cm2 + offset.g.cm2
max.depth.g.cm2         <- max.depth.g.cm2 + offset.g.cm2

# read in DDD files  
ddds                    <- dataDDDset(ddd.path = ddd.path)

# from those, select IESs necessary
jj                      <- ddds@peak.positions.g.cm2 > min.depth.g.cm2 & ddds@peak.positions.g.cm2 <= max.depth.g.cm2 &
                             rep(c(TRUE, rep(FALSE, IES.step - 1)), length.out = length(ddds@projectiles))
jj[head(which(jj), extend.IES[1]) - IES.step] <- TRUE
jj[tail(which(jj), extend.IES[2]) + IES.step] <- TRUE
no.IES                          <- sum(jj)
ddds.sub                        <- ddds[which(jj)]


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

  # Get normalized deviation
  dev <- (doses - dose.set.Gy)/dose.set.Gy
  return(sum(dev^2))
}

# minimize objective function to find weights
if(.Platform$OS.type == "unix") {
  cl                <- makeForkCluster()
  rel.weights       <- optimParallel( fn             = dose.dev, 
                                      par            = rep( 1, no.IES) + seq(0, 1, length.out = no.IES),
                                      depths.g.cm2   = depths.g.cm2,
                                      DDD.set        = ddds.sub,
                                      dose.set.Gy    = plateau.dose.per.primary.Gy,
                                      method         = "L-BFGS-B",
                                      lower          = rep(0.1, no.IES),
                                      upper          = rep(20, no.IES),
                                      control        = list(trace = 1, 
                                                            maxit = 300,
                                                            factr = conversion.threshold/.Machine$double.eps),
                                      parallel       = list(cl = cl, forward = FALSE, loginfo = FALSE))$par
  stopCluster(cl)
} else {
  rel.weights       <- optim        ( fn             = dose.dev, 
                                      par            = rep( 1, no.IES) + seq(0, 1, length.out = no.IES),
                                      depths.g.cm2   = depths.g.cm2,
                                      DDD.set        = ddds.sub,
                                      dose.set.Gy    = plateau.dose.per.primary.Gy,
                                      method         = "L-BFGS-B",
                                      lower          = rep(0.1, no.IES),
                                      upper          = rep(20, no.IES),
                                      control        = list(trace = 1, 
                                                            maxit = 300,
                                                            factr = conversion.threshold/.Machine$double.eps))$par
  
  
}
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


####################################
# B. Read spc files needed if chosen
if((biol.optimization & (rbe.method == "SPCs and RBE file"))| output.spectra){
  #################
  # 1. Read in data
  
  # Scan available spc files to get energy sample points
  spcs             <- dataSPCset(spc.path = spc.path)
  
  # Returns a list of dataSPC objects with same energies
  # as selected IESs. To this end, dataSPCs will be interpolated
  # between the adjacent energies of the available
  # energy sample point in the spc set
  ddds.spc          <- mclapply(ddds.sub@beam.energies.MeV.u, 
                                function(x, s){
                                  get.spc(s, x)},
                                spcs)

  rbe.data <- dataRBE(rbe.file, rbe.path)
}  

############################
# C. Biological optimization
if(biol.optimization & (rbe.method == "SPCs and RBE file")){
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
                                control        = list(trace = 1, 
                                                      maxit = 300,
                                                      factr = conversion.threshold/.Machine$double.eps))$par
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

}

if(biol.optimization & (rbe.method == "From alpha/beta with depth")){
  ###############################
  # Iterate with RBE weighting
  for(biol.opt.step in 1:n.biol.opt.steps){
    # biol.opt.step <- 1
    
    D.phys.Gy             <- get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                   depths.g.cm2 = depths.g.cm2, 
                                                   weights      = total.weights)
    rbe                   <- get.RBE.from.set(     DDD.set      = ddds.sub, 
                                                   depths.g.cm2 = depths.g.cm2,
                                                   weights      = total.weights)
    
    plateau.dose.per.primary.Gy <- mean(get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                                             depths.g.cm2 = depths.g.cm2))
    # Scale down using RBE at all depth of SOBP
    rbe.scaled.plateau.doses.per.primary.Gy <- plateau.dose.per.primary.Gy / rbe
    fluence.factor <- plateau.dose.Gy / mean(plateau.dose.per.primary.Gy)
    
    # minimize objective function to find weights, use current weights as starting parameters
    rel.weights         <- optim( fn             = dose.dev, 
                                  par            = total.weights / fluence.factor,
                                  depths.g.cm2   = depths.g.cm2,
                                  DDD.set        = ddds.sub,
                                  dose.set.Gy    = rbe.scaled.plateau.doses.per.primary.Gy,
                                  method         = "L-BFGS-B",
                                  lower          = rep(0.0, no.IES),
                                  control        = list(trace = 1, 
                                                        maxit = 300,
                                                        factr = conversion.threshold/.Machine$double.eps))$par
    # Scale number of primaries to get actual weights
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
  
  D.phys.Gy             <- get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                 depths.g.cm2 = bio.depths.g.cm2, 
                                                 weights      = total.weights)
  
  rbe                   <- get.RBE.from.set(     DDD.set      = ddds.sub, 
                                                 depths.g.cm2 = bio.depths.g.cm2,
                                                 weights      = total.weights)
  
  df.plot <- data.frame( depth.g.cm2          = plot.depths.g.cm2,
                         D.phys.Gy            = get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                                      depths.g.cm2 = plot.depths.g.cm2, 
                                                                      weights      = total.weights),
                         RBE                  = approx( x    = bio.depths.g.cm2,
                                                        y    = rbe,
                                                        xout = plot.depths.g.cm2,
                                                        rule = 2)$y)
  df.plot$D.biol.Gy <- df.plot$D.phys.Gy * df.plot$RBE

}


###################
# D. Write SOBP.dat
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

if(biol.optimization){
  myyscale.component <- function(...)
  {
    ans 				              <- yscale.components.default(...)
    ans$right 	              <- ans$left
    foo 				              <- ans$right$labels$at
    ans$right$labels$labels 	<- as.character(foo * 5)
    ans
  }
  
  xyplot(D.biol.Gy + D.phys.Gy + RBE/5 ~ depth.g.cm2,
         df.plot,
         type = "l",
         xlab = list("depth / (g/cm2)", cex=1.5),
         ylab = list("Dose / (Gy, GyRBE)", cex=1.5),
         par.settings=	list( layout.widths=list(right.padding=10)),
         yscale.component	= myyscale.component,
         legend = list(right = list(fun = textGrob,
                                    args = list(x = 3,
                                                y = 0.5,
                                                rot = 90,
                                                label = "RBE",
                                                gp=gpar(cex = 1.5, col = "darkgreen"),
                                                just = "center",
                                                default.units = "native",
                                                vp = viewport(xscale = c(0, 1),
                                                              yscale = c(0, 1))))),
         scales = list( x = list(cex = 1.25),
                        y = list(cex = 1.25,
                                 relation = "free",
                                 rot = 0)),
         grid = TRUE,
         main = list(paste0("SOBP (single field, ",
                            unique(ddds.sub@projectiles),
                            ") consisting of ",
                            no.IES,
                            " IESs "),
                     cex=1.5),
         sub  = "NB: HIT isocenter is at 0.289 g/cm2 depth!")
}

################################
# E. Write out spectra, do plots
if(output.spectra & (rbe.method == "SPCs and RBE file")){
  # Generate vector with depth positions
  spectra.depths.g.cm2 <- seq(0, max.depth.g.cm2*2, by = spectra.step.size.g.cm2)
  
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
                                d = spectra.depths.g.cm2)
  
  # Weights spectra by fluence weights of SOBP
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
  
  # Adds weighted spectra
  eff.spectra.at.depth <- mclapply(1:length(spectra.depths.g.cm2),
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
  
  # Create spectra data table for output
  dt.spectra <- get.data.table(eff.spectra.at.depth)
  write.table(dt.spectra, file = paste0(expid, "_spectra.csv"), quote = FALSE, sep = "; ", row.names = FALSE)
}

alarm()
