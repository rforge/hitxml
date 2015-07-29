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
ddd.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/12C/RF3MM"

# minimal and maximal depth in cm
min.depth         <- 10
max.depth         <- 15


#############################
# read in DDD files  

ddds             <- dataDDDset(ddd.path = ddd.path)

#############################
# calculate the number of 
# available IESs and equalize 
# their maximum water depth
jj                <- which(ddds@peak.positions.g.cm2 > min.depth & ddds@peak.positions.g.cm2 < max.depth)

no.IES            <- length(jj)
ddds.sub          <- ddds[jj]

z.max             <- min.depth + max.depth


##############
# model SOBP #
##############

resolution         <- 100
plateau.dose.Gy    <- 2

z.max              <- round(z.max, digits=1)
depth.seq          <- seq(from=0, to=z.max, by=z.max/resolution)

min.depth.step     <- ceiling( head(sort(ddds.sub@peak.positions.g.cm2), 1) / (z.max/resolution) )
max.depth.step     <- floor( tail(sort(ddds.sub@peak.positions.g.cm2), 1) / (z.max/resolution) )

get.dose.Gy.from.set(ddds.sub, 5)


primary.weights    <- abs(nlm(f=HX.deviation.proton, p=c(1:no.IES))$estimate)
total.weights      <- round( primary.weights * plateau.dose.Gy / mean(HX.SOBP.dose(primary.weights)[min.depth.step:max.depth.step]), digits=0 )
used.beam.energies <- beam.energies[jj]

##############
# save SOBP  #
##############
n.data             <- length(used.beam.energies)
df                 <- data.frame(E.GeV   = used.beam.energies / 1000,
                                 x.cm    = rep(0.0, n.data),
                                 y.cm    = rep(0.0, n.data),
                                 FWHM.cm = rep(0.1, n.data),
                                 fluence = as.integer(total.weights))
write.table(df,
            file = "SOBP.dat",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol       = "\r\n" )

############################
# plot SOBP (single field) #
############################

plot1  <- xyplot(HX.SOBP.dose(total.weights) ~ depth.seq,
                 type = "l",
                 lwd  = 3,
                 col ="blue",
                 #xlim = c(0, z.max),
                 #ylim = c(0, 1.1*max(SOBP.dose(total.weights))),
                 xlab = list("distal depth [cm]", cex=1.5),
                 ylab = list("total dose [Gy]", cex=1.5),
                 main = list(paste("SOBP (single field)", "consisting of", no.IES, "IESs"), cex=1.5))

# layer1 <- layer( panel.segments(x0=min.depth.cm, x1=min.depth.cm, y0=0, y1=1.1*max(SOBP.dose(total.weights)), lwd=2, col="black"),
#                  panel.segments(x0=max.depth.cm, x1=max.depth.cm, y0=0, y1=1.1*max(SOBP.dose(total.weights)), lwd=2, col="black") )

plot(plot1)