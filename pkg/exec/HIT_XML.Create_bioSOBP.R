############################
# Script to generate an
# biologically opt. SOBP
############################
# Original version by
# Grischa Klimpki, 2014
# adapted to HIT_XML
# package by Steffen
# Greilich, 2015-07
############################
rm(list = ls())

library(HITXML)

###########################
# define input parameters #

# path to spc and ddd data
ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/12C/RF3MM"
spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM"
rbe.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/RBE"

# minimal and maximal depth in cm
min.depth.g.cm2         <- 10
max.depth.g.cm2         <- 15

step.size.g.cm2         <- 0.1
IES.step                <- 3
plateau.dose.Gy         <- 2.0

biol.optimization       <- TRUE
rbe.file                <- "target02.rbe"
n.biol.opt.steps        <- 5
bio.step.size.g.cm2     <- 0.25

write.SOBP              <- FALSE

plot.range              <- 1.5

#############################
# read in DDD files  

ddds             <- dataDDDset(ddd.path = ddd.path)


#############################
# calculate the number of 
# available IESs and equalize 
# their maximum water depth
jj                <- ddds@peak.positions.g.cm2 >= min.depth.g.cm2 & ddds@peak.positions.g.cm2 <= max.depth.g.cm2 &
                             rep(c(TRUE, rep(FALSE, IES.step - 1)), length.out = length(ddds@projectiles))
jj[head(which(jj), 1)-IES.step] <- TRUE
jj[tail(which(jj), 1)+IES.step] <- TRUE
no.IES            <- sum(jj)

ddds.sub          <- ddds[which(jj)]

###### model physical SOBP

# Vector of depths covering the SOBP
depths.g.cm2       <- seq( from       = min.depth.g.cm2, 
                           to         = max.depth.g.cm2, 
                           by         = step.size.g.cm2)


# Objective function: sum of squares of deviation to dose set
mean.dose.Gy <- mean(get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                          depths.g.cm2 = depths.g.cm2))

dose.dev <- function(p, depths.g.cm2, DDD.set, dose.set.Gy){
  doses <- get.dose.Gy.from.set(DDD.set      = DDD.set, 
                                depths.g.cm2 = depths.g.cm2, 
                                weights      = p)
  log10(sum( (doses - dose.set.Gy)^10))
}


# minimize objective function to find weights
rel.weights       <- optim( fn             = dose.dev, 
                            par            = rep( 1, no.IES),
                            depths.g.cm2   = depths.g.cm2,
                            DDD.set        = ddds.sub,
                            dose.set.Gy    = mean.dose.Gy,
                            method         = "L-BFGS-B",
                            lower          = rep(0.0, no.IES),
#                            upper          = rep(20, no.IES),
                            control        = list(trace = TRUE, 
                                                  maxit = 200))$par

total.weights     <- rel.weights * plateau.dose.Gy / mean.dose.Gy / no.IES

##########################
# plot SOBP (single field) 
pt.depths.g.cm2         <- seq( from       = 0.0, 
                                to         = max(ddds.sub@peak.positions.g.cm2) * plot.range, 
                                by         = step.size.g.cm2)

plot.SOBP(ddds.sub, 
          pt.depths.g.cm2, 
          total.weights, 
          "(phys.opt.)")


#################################
# B. Biological optimzation steps
if(biol.optimization){
  #################
  # 1. Read in data
  
  # Scan available spc files to get energy sample points
  spcs             <- dataSPCset(spc.path = spc.path)
  
  # Returns a list of dataSPC objects with same energies
  # than selected IESs. dataSPCs will be interpolated
  # between the adjacent energies of the available
  # energy sample point in the spc set
  ddds.spc          <- lapply(ddds.sub@beam.energies.MeV.u, 
                              function(x, s){
                                get.spc(s, x)},
                              spcs)
  # save(ddds.spc, file = "ddds.spc.rda")
  # load("ddds.spc.rda")
  
  # no.IES <- 1
  # depths.g.cm2 <- 0
  
  # Returns a list (length number of IESs in SOBP)
  # of lists (each length number of depths covering SOBP)
  # of dataSpectrum objects - representing the spectra
  # from the individual IESs (first index) contributing at specific depth (second index)
  spectra.at.depth  <- lapply(1:no.IES,
                              function(i, s, d){
                                cat("Getting spectra from IES", i, "\n")
                                dataSpectrum(s[[i]], d)},
                              s = ddds.spc,
                              d = depths.g.cm2)
  # save(spectra.at.depth, file = "sad.rda")
  # load("sad.rda")

  rbe.data <- dataRBE(rbe.file, rbe.path)
  
  ###############################
  # 2. Iterate with RBE weighting
  for(biol.opt.step in 1:n.biol.opt.steps){
    # Applies the given weights to all spectra at depths covering SOBP (second index)
    # from respective IESs (first index) 
    w.spectra.at.depth  <- lapply(1:no.IES,
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
    eff.spectra.at.depth <- lapply(1:length(depths.g.cm2),
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
                                                   weights      = total.weights)*no.IES
    rbe      <- HX.RBE.LEM(rbe.data, eff.spectra.at.depth, D.phys.Gy)$RBE
    
    ## SECOND RUN: replace by functions!
    mean.dose.Gy <- mean(get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                              depths.g.cm2 = depths.g.cm2)) / rbe
    
    # minimize objective function to find weights
    rel.weights       <- optim( fn             = dose.dev, 
                                par            = rep( 1, no.IES),
                                depths.g.cm2   = depths.g.cm2,
                                DDD.set        = ddds.sub,
                                dose.set.Gy    = mean.dose.Gy,
                                method         = "L-BFGS-B",
                                lower          = rep(0.0, no.IES),
                                control        = list(trace = TRUE, 
                                                      maxit = 200))$par
    
    total.weights     <- rel.weights * plateau.dose.Gy / mean(mean.dose.Gy * rbe)  / no.IES
    
    plot(plot.SOBP(ddds.sub, 
                   pt.depths.g.cm2, 
                   total.weights, 
                   paste0("(biol.opt., step ", biol.opt.step, ")")))
  }

  ###########
  # GET RBE for entire field ( with reduced resolution)
  bio.depths.g.cm2         <- seq( from       = 0.0, 
                                  to         = max(ddds.sub@peak.positions.g.cm2) * plot.range, 
                                  by         = bio.step.size.g.cm2)
  
  spectra.at.depth  <- lapply(1:no.IES,
                              function(i, s, d){
                                cat("Getting spectra from IES", i, "\n")
                                dataSpectrum(s[[i]], d)},
                              s = ddds.spc,
                              d = bio.depths.g.cm2)

  w.spectra.at.depth  <- lapply(1:no.IES,
                                function(i, s, w){
                                  cat("Weighting spectra from IES", i, "\n")
                                  lapply(1:length(s[[i]]),
                                         function(j, ss, ww){
                                           ss[[j]] * ww},
                                         ss = s[[i]],
                                         ww = w[i])},
                                s = spectra.at.depth,
                                w = total.weights)
  
  eff.spectra.at.depth <- lapply(1:length(bio.depths.g.cm2),
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
                                                 weights      = total.weights)*no.IES
  
  rbe      <- HX.RBE.LEM(rbe.data, eff.spectra.at.depth, D.phys.Gy)$RBE

  df.plot <- data.frame( depth.g.cm2          = pt.depths.g.cm2,
                         D.phys.Gy            = get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                 depths.g.cm2 = pt.depths.g.cm2, 
                                                 weights      = total.weights)*no.IES,
                         RBE                  = approx( x    = bio.depths.g.cm2,
                                                        y    = rbe,
                                                        xout = pt.depths.g.cm2,
                                                        rule = 2)$y)
  df.plot$D.biol.Gy <- df.plot$D.phys.Gy * df.plot$RBE
  
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
                                                gp=gpar(cex = 1.5),
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
                     cex=1.5))
}

###############
# save SOBP.dat
if(write.SOBP){
  df                 <- data.frame(E.GeV   = ddds.sub@beam.energies.MeV.u / 1000,
                                   x.cm    = rep(0.0, no.IES),
                                   y.cm    = rep(0.0, no.IES),
                                   FWHM.cm = rep(0.1, no.IES),
                                   fluence = as.integer(total.weights))
  write.table(df,
              file = "SOBP.dat",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              eol       = "\r\n" )
}

