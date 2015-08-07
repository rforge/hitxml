rm(list = ls())
library(HITXML)
library(libamtrack)
library(lattice)

# TO SET: meta information
projectile            <- "12C"
beam.energy.MeV.u     <- 270.00
target.material       <- "Water, Liquid"
peak.position.g.cm2   <- 14.08
ripple.filter         <- c("0MM", "3MM")[2]

# TO SET: number of primaries
n.primaries           <- 5 * 5 * 100

# TO SET: maximum energy and spacing (linear)
min.E.MeV.u           <- 0.0
max.E.MeV.u           <- 600.0
n.bins                <- 600
log.bins              <- FALSE

# TO SET: material and stopping power data
material.name         <- "Water, Liquid"
stopping.power.source <- "ICRU"

# TO SET: data location
data.path             <- "."



# READ FLUKA OUTPUT (can be multiple files)
df <- read.FLUKA.phase.space(data.path, Z.max = 6)
save(df, file = "df.rda")
# load("df.rda")

depth.steps   <- 1:79
depths.g.cm2  <- c(0.0000, 0.0500, 0.1000, 0.1500, 0.2000, 
                   0.2500, 0.3000, 0.3500, 0.4000, 0.4500, 
                   0.5000, 0.5500, 0.6000, 0.6500, 0.7000, 
                   0.7500, 0.8000, 0.8100, 0.8200, 0.8300, 
                   0.8400, 0.8500, 0.8600, 0.8700, 0.8800, 
                   0.8900, 0.9000, 0.9100, 0.9200, 0.9300, 
                   0.9400, 0.9500, 0.9550, 0.9600, 0.9650, 
                   0.9700, 0.9750, 0.9800, 0.9850, 0.9900, 
                   0.9910, 0.9920, 0.9930, 0.9940, 0.9950, 
                   0.9960, 0.9970, 0.9980, 0.9990, 0.9995, 
                   1.0000, 1.0005, 1.0010, 1.0020, 1.0030, 
                   1.0040, 1.0050, 1.0060, 1.0070, 1.0080, 
                   1.0090, 1.0100, 1.0120, 1.0140, 1.0160, 
                   1.0180, 1.0200, 1.0250, 1.0300, 1.0400, 
                   1.0500, 1.0600, 1.0700, 1.0800, 1.0900, 
                   1.1000, 1.1100, 1.4100, 1.9100) * peak.position.g.cm2

# Get boundaries
cv        <- unique(df$travel)

# Convert
ee        <- lapply(cv,
                    function(cur.cv){ cat(paste0("Batch: ", cur.cv, "\n"))
                      ii <- df$travel == cur.cv
                      hh <- get.multivariate.weighted.histogram( df$EKIN[ii],
                                                                 df$ICHRGE[ii],
                                                                 1/df$CZTRCK[ii],
                                                                 min.E.MeV.u,
                                                                 max.E.MeV.u,
                                                                 n.bins,
                                                                 log.bins)
                      names(hh) <- c("Z", "frequency", "E.MeV.u", "dE.MeV.u", "density")
                      hh$travel <- cur.cv
                      hh})

# Get particles per primary
ee             <- lapply(1:length(ee),
                         function(idx){ tmp    <- as.data.frame(ee[[idx]])
                                        tmp$frequency / n.primaries
                                        tmp$density   / n.primaries
                                        tmp})

dfE <- do.call("rbind.data.frame", ee)

# Add boundary depth
dfE$depth.step  <- depth.steps[match(dfE$travel, cv)]
dfE$depth.g.cm2 <- depths.g.cm2[match(dfE$travel, cv)]

save(dfE, file = "dfE.sg74000.Rda")
# load("dfE.sg74000.Rda")

# Create SPC
spc <- dataSPC(projectile          = projectile,
               beam.energy.MeV.u   = beam.energy.MeV.u, 
               target.material     = target.material, 
               peak.position.g.cm2 = peak.position.g.cm2, 
               depth.step          = dfE$depth.step,
               depth.g.cm2         = dfE$depth.g.cm2, 
               particle.no         = dfE$Z * 1002, 
               E.MeV.u             = dfE$E.MeV.u, 
               dE.MeV.u            = dfE$dE.MeV.u, 
               N.per.primary       = dfE$frequency)

name.spc     <- paste(projectile, "generic", beam.energy.MeV.u, ripple.filter, sep = "_")
save(spc, file = paste0(name.spc, ".Rda"))

# FINISHED, plot results
plot.comment <- paste0( projectile, 
                        " on ", 
                        target.material, 
                        " (", 
                        beam.energy.MeV.u, 
                        " MeV/u), ripple filter: ",
                        ripple.filter)

xyplot(frequency ~ E.MeV.u|sprintf("Depth %03d mm", depth.g.cm2 * 10),
       dfE,
       type   = "s",
       groups = Z,
       auto.key = list(space = "right", lines = TRUE, points = FALSE),
       scales = list( y = list(log = 10)),
       grid   = TRUE,
       as.table = TRUE,
       xlab     = "kinetic energy / (MeV/amu)",
       main     = plot.comment,
       layout   = c(5,4))

