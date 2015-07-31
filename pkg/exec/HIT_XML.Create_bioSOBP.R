############################
# Script to generate an
# SOBP
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
#ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/p/RF0MM"
ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/12C/RF3MM"
spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM"

# minimal and maximal depth in cm
min.depth.g.cm2         <- 10
max.depth.g.cm2         <- 15

step.size.g.cm2         <- 0.1
IES.step                <- 3
plateau.dose.Gy         <- 2.0

write.SOBP              <- FALSE

#############################
# read in DDD files  

ddds             <- dataDDDset(ddd.path = ddd.path)
spcs             <- dataSPCset(spc.path = spc.path)


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

ddds.spc          <- lapply(ddds.sub@beam.energies.MeV.u,
                            function(x, s){
                              get.spc(s, x)},
                            spcs)

# save(ddds.spc, file = "ddds.spc.rda")
# load("ddds.spc.rda")
############
# model SOBP

depths.g.cm2       <- seq( from       = min.depth.g.cm2, 
                           to         = max.depth.g.cm2, 
                           by         = step.size.g.cm2)

weights           <- 1:no.IES
spectra.at.depth  <- lapply(1:no.IES,
                            function(i, s, d){
                            cat("Getting spectra from IES", i, "\n")
                            cur.spectra <- dataSpectrum(s[[i]], d)},
                            s = ddds.spc,
                            d = depths.g.cm2)
# save(spectra.at.depth, file = "sad.rda")
w.spectra.at.depth  <- lapply(1:no.IES,
                              function(i, s, w){
                                cat("Weighting spectra from IES", i, "\n")
                                lapply(1:length(s[[i]]),
                                       function(j, ss, ww){
                                         ss[[j]] * ww[i]},
                                       ss = s[[i]],
                                       ww = w)},
                              s = spectra.at.depth,
                              w = weights)

pp <- lapply(1:length(depths.g.cm2),
             function(i, s, n, w){
               cat("Adding spectra for depth", i, "\n")
               ss <- s[[1]][[i]]
               for(j in 2:n){
                 ss <- ss + s[[j]][[i]]
               }
               ss},
               s = spectra.at.depth,
               n = no.IES)

               plot(pp[[50]])
               dose.from.spectrum.Gy(pp[[50]], stopping.power.source = "ICRU", target.material = "Water, Liquid")
               
               
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

##########################
# plot SOBP (single field) 
depths.g.cm2       <- seq( from       = 0.0, 
                           to         = max(ddds.sub@peak.positions.g.cm2) * 1.75, 
                           by         = step.size.g.cm2)

D.Gy               <- get.dose.Gy.from.set(DDD.set      = ddds.sub, 
                                           depths.g.cm2 = depths.g.cm2, 
                                           weights      = total.weights)*no.IES
xyplot(D.Gy ~ depths.g.cm2,
       grid = TRUE,
       type = "l",
       xlab = list("distal depth [cm]", cex=1.5),
       ylab = list("total dose [Gy]", cex=1.5),
       main = list(paste("SOBP (single field)", "consisting of", no.IES, "IESs"), cex=1.5))

dd <- unlist(lapply(1:no.IES,
       function(i){
         get.dose.Gy(get.ddd(ddds.sub, ddds.sub@beam.energies.MeV.u[i]), depths.g.cm2) * total.weights[i]
       }))

df <- data.frame( depths.g.cm2 = rep(depths.g.cm2, no.IES+1),
                  D.Gy         = c(D.Gy, dd),
                  which        = c(rep("SOBP", length(depths.g.cm2)),
                                   sort(rep(1:no.IES, length(depths.g.cm2)))))
xyplot(D.Gy ~ depths.g.cm2,
       df,
       grid = TRUE,
       type = "l",
       col  = c(rep("grey", no.IES), "blue"),
       alpha  = c(rep(0.8, no.IES), 1),
       groups = which,
       xlab = list("distal depth [cm]", cex=1.5),
       ylab = list("total dose [Gy]", cex=1.5),
       main = list(paste0("SOBP (single field, ",
                          unique(ddds.sub@projectiles),
                          ") consisting of ", no.IES, " IESs"), 
                   cex=1.5))

