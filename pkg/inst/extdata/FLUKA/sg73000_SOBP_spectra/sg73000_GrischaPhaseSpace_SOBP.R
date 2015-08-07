rm(list = ls())
library(HITXML)
library(libamtrack)
library(lattice)

# TO SET: meta information
plot.comment          <- "sg73000 (SOBP C-12, 10-15 cm)"
boundary.depth.cm     <- c(1,8,11,13.1,14,17)
data.path             <- "./data_sg73000_SOBP"

# TO SET: maximum LET (keV/um) and spacing (linear)
min.LET.keV.um        <- 0.1
max.LET.keV.um        <- 300.0
n.bins                <- 300
log.bins              <- TRUE

# TO SET: material and stopping power data
material.name         <- "Water, Liquid"
stopping.power.source <- "ICRU"

# TO SET: dose at boundary with index idx
dose.Gy               <- 2.0
dose.idx              <- 4


# READ FLUKA OUTPUT (can be multiple files)
df <- read.FLUKA.phase.space(data.path, Z.max = 6)
# save(df, file = "df.rda")
# load("df.rda")

# Get boundaries
cv        <- unique(df$travel)

# Convert, all dose = 1
ll        <- lapply(cv,
                    function(cur.cv){ cat(paste0("Batch: ", cur.cv, "\n"))
                                 ii <- df$travel == cur.cv
                                 hh <- get.multivariate.weighted.histogram( df$LET.UNR[ii],
                                                                            df$ICHRGE[ii],
                                                                            1/df$CZTRCK[ii],
                                                                            min.LET.keV.um,
                                                                            max.LET.keV.um,
                                                                            n.bins,
                                                                            TRUE)
                                 names(hh) <- c("Z", "frequency", "LET.H2O.keV.um", "dLET", "density")
                                 hh$travel <- cur.cv
                                 hh})

ee        <- lapply(cv,
                    function(cur.cv){ cat(paste0("Batch: ", cur.cv, "\n"))
                      ii <- df$travel == cur.cv
                      hh <- get.multivariate.weighted.histogram( df$EKIN[ii],
                                                                 df$ICHRGE[ii],
                                                                 1/df$CZTRCK[ii],
                                                                 0,
                                                                 300,
                                                                 n.bins,
                                                                 FALSE)
                      names(hh) <- c("Z", "frequency", "E.MeV.u", "dE.MeV.u", "density")
                      hh$travel <- cur.cv
                      hh})

# Get relative doses
doses          <- sapply(ll, function(x){sum(x$LET.H2O.keV.um * x$frequency * 1.602e-10)})
fluence.factor <- (doses / doses[dose.idx]) * (dose.Gy / doses[dose.idx])

ll             <- lapply(1:length(ll),
                         function(idx){ tmp    <- as.data.frame(ll[[idx]])
                                        tmp$frequency * fluence.factor[idx]
                                        tmp$density   * fluence.factor[idx]
                                        tmp})

dfL <- do.call("rbind.data.frame", ll)
dfE <- do.call("rbind.data.frame", ee)

dfL$Z <- as.character(dfL$Z)
cvT  <- paste(dfL$LET.H2O.keV.um, dfL$travel)
dfL   <- rbind(dfL,
               data.frame( Z              = tapply(dfL$Z,              cvT, function(x){"total"}),
                           frequency      = tapply(dfL$frequency,      cvT, sum),
                           LET.H2O.keV.um = tapply(dfL$LET.H2O.keV.um, cvT, unique),
                           dLET           = tapply(dfL$dLET,           cvT, unique),
                           density        = tapply(dfL$density,        cvT, sum),
                           travel         = tapply(dfL$travel,         cvT, unique)))
                           
dfL   <- dfL[order(dfL$travel, dfL$Z, dfL$LET.H2O.keV.um),]

# Get depth of boundary
dfL$depth.cm <- boundary.depth.cm[match(dfL$travel, cv)]
dfE$depth.cm <- boundary.depth.cm[match(dfE$travel, cv)]

save(dfL, file = "dfL.sg73000.Rda")
save(dfE, file = "dfE.sg73000.Rda")

# FINISHED, plot results
xyplot(frequency ~ LET.H2O.keV.um|sprintf("Depth %03d mm", depth.cm * 10),
       dfL,
       type   = "l",
       groups = Z,
       auto.key = list(space = "right", lines = TRUE, points = FALSE),
       scales = list(x = list(log = 10), y = list(relation = "free")),
       grid   = TRUE,
       as.table = TRUE,
       xlab     = "LET in water / (MeV.cm)",
       main     = plot.comment)

xyplot(frequency ~ E.MeV.u|sprintf("Depth %03d mm", depth.cm * 10),
       dfE,
       type   = "s",
       groups = Z,
       auto.key = list(space = "right", lines = TRUE, points = FALSE),
       scales = list( y = list(log = 10)),
       grid   = TRUE,
       as.table = TRUE,
       xlab     = "kinetic energy / (MeV/amu)",
       main     = plot.comment)

