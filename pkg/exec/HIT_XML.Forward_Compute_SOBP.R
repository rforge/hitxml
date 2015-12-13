#' Script to forward compute a spread-out Bragg peak
#'
#' Steffen Greilich, 2015-12

rm(list = ls())
library(HITXML)

#' START OF USER INPUT
SOBP.path <- "C:/Users/greilich/Desktop/SOBP.dat"
  
# path to spc and ddd data
ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/12C/RF3MM"
spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM"
rbe.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/RBE"

depth.min.g.cm2          <- 0.0
depth.max.g.cm2          <- 20.0
step.size.g.cm2          <- 0.05

biol.dose                <- TRUE
rbe.file                 <- "target02.rbe"
bio.step.size.g.cm2      <- 0.25
spectra.at.depths.g.cm2  <- c(1,2,3)

#' END OF USER INPUT

# read in SOBP file
df.SOBP   <- read.table(SOBP.path,
                        header = FALSE,
                        col.names = c("E.GeV", "x.cm", "y.cm", "FWHM.cm", "fluence"))

# read in DDD files  
ddds              <- dataDDDset(ddd.path = ddd.path)

# get IESs necessary
jj                <- round(ddds@beam.energies.MeV.u, 2)%in%round(df.SOBP$E.GeV*1000, 2)
ddds.sub          <- ddds[which(jj)]
if(length(ddds.sub@beam.energies.MeV.u) != nrow(df.SOBP)){
  stop("Something funny - number of IESs does not match...")
}

#' A. DO PHYSICAL DOSE

# Vector of depths covering the SOBP
depths.g.cm2       <- seq( from       = depth.min.g.cm2, 
                           to         = depth.max.g.cm2, 
                           by         = step.size.g.cm2)
bio.depths.g.cm2       <- seq( from       = depth.min.g.cm2, 
                                to         = depth.max.g.cm2, 
                                by         = bio.step.size.g.cm2)

doses <- get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                              depths.g.cm2 = depths.g.cm2, 
                              weights      = )

plot.SOBP(plot.ddds = ddds.sub, 
          plot.depths.g.cm2 = depths.g.cm2, 
          plot.weights = df.SOBP$fluence, 
          "(phys.opt.)")
  

# B. Biological dose
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
  # TODO: Move this format to its own class and adapt lapply's below
  spectra.at.depth  <- lapply(1:nrow(df.SOBP),
                              function(i, s, d){
                                cat("Getting spectra from IES", i, "\n")
                                dataSpectrum(s[[i]], d)},
                              s = ddds.spc,
                              d = bio.depths.g.cm2)
  # save(spectra.at.depth, file = "sad.rda")
  # load("sad.rda")

  rbe.data <- dataRBE(rbe.file, rbe.path)
  
 # Applies the given weights to all spectra at depths covering SOBP (second index)
    # from respective IESs (first index) 
    w.spectra.at.depth  <- lapply(1:nrow(df.SOBP),
                                  function(i, s, w){
                                    cat("Weighting spectra from IES", i, "\n")
                                    lapply(1:length(s[[i]]),
                                           function(j, ss, ww){
                                             ss[[j]] * ww},
                                           ss = s[[i]],
                                           ww = w[i])},
                                  s = spectra.at.depth,
                                  w = df.SOBP$fluence)
    
    # Combines (adds) spectra from all IESs (first index) for the depths 
    # covering the SOBP (second index). Spectra can be have been weighted before
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
                 n = nrow(df.SOBP))
    
    # Why does this get different results??
    #dose.Gy(pp, 
    #        stopping.power.source = "ICRU", 
    #        target.material = "Water, Liquid")
    D.phys.Gy             <- get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                   depths.g.cm2 = bio.depths.g.cm2, 
                                                   weights      = df.SOBP$fluence)
    rbe                   <- HX.RBE.LEM(rbe.data, 
                                        eff.spectra.at.depth, 
                                        D.phys.Gy)$RBE
    

  df.plot <- data.frame( depth.g.cm2          = depths.g.cm2,
                         D.phys.Gy            = get.dose.Gy.from.set( DDD.set      = ddds.sub, 
                                                                      depths.g.cm2 = depths.g.cm2, 
                                                                      weights      = df.SOBP$fluence),
                         RBE                  = approx( x    = bio.depths.g.cm2,
                                                        y    = rbe,
                                                        xout = depths.g.cm2,
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
                            nrow(df.SOBP), 
                            " IESs "), 
                     cex=1.5),
         sub  = "NB: HIT isocenter is at 0.289 g/cm2 depth!")
}



#############
#

depth <- lapply(eff.spectra.at.depth,
                spectrum.get.depth.g.cm2)
total.fluence <- lapply(eff.spectra.at.depth,
                        spectrum.total.n.particles)
dose.Gy <-lapply(eff.spectra.at.depth,
                 spectrum.dose.Gy,
                 "ICRU",
                 "Water, Liquid")
fLET <-lapply(eff.spectra.at.depth,
              function(x){
                LET <- spectrum.Mass.Stopping.Power.MeV.cm2.g(x,
                 "ICRU",
                 "Water, Liquid")
                return(sum(x@spectrum[,"N"]*LET)/sum(x@spectrum[,"N"]))}
              )
dLET <-lapply(eff.spectra.at.depth,
              function(x){
                LET <- spectrum.Mass.Stopping.Power.MeV.cm2.g(x,
                                                              "ICRU",
                                                              "Water, Liquid")
                return(sum(x@spectrum[,"N"]*LET*LET)/sum(x@spectrum[,"N"]*LET))}
)

xyplot(total.fluence ~ depth)
xyplot(dose.Gy ~ depth)
xyplot(fLET ~ depth)
xyplot(dLET ~ depth)

spectra.idx <- match(spectra.at.depths.g.cm2, depth)

sss <- eff.spectra.at.depth[spectra.idx]
xyplot(sss[[3]]@spectrum[,"N"]~sss[[3]]@spectrum[,"E.MeV.u"],
       groups = sss[[3]]@spectrum[,"particle.no"],
       auto.key = TRUE,
       scales = list(y = list(log = 10)))