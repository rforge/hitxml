############################
# Script to generate an
# SOBP
############################
# Original version by
# Grischa Klimpki, 2014
# adapted to HIT_XML
# package by Steffen
# Greilich, 2015-07-15
############################
rm(list = ls())

library(libamtrack)
library(lattice)

source("./exec/Generate_SOBP_Functions.R")

###########################
# define input parameters #
###########################

# path to spc and ddd data
data.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/12C/RF3MM"

# minimal and maximal depth in cm
min.depth         <- 10
max.depth         <- 15


#############################
# read in DDD files  
#############################

DDD.file.names   <- list.files(path = data.path, full.names=TRUE)
split            <- strsplit(DDD.file.names, paste0(data.path, "/12C_E"))
split            <- as.numeric(sapply(split, function(x) x <- sub("_rifi3mm_ddd_new.ddd", "", x[2])))
DDD.file.names   <- DDD.file.names[order(split)]

beam.energies    <- vector(mode="list", length=length(DDD.file.names))
beam.energies    <- as.numeric(lapply(1:length(DDD.file.names), read.beam.energy, file.names=DDD.file.names))

BP.depth         <- vector(mode="list", length=length(DDD.file.names))
BP.depth         <- as.numeric(lapply(1:length(DDD.file.names), get.BP.depth, file.names=DDD.file.names))

DDD.list         <- vector(mode="list", length=length(DDD.file.names))
DDD.list         <- lapply(1:length(DDD.file.names), read.DDD, file.names=DDD.file.names)

#############################
# calculate the number of 
# available IESs and equalize 
# their maximum water depth
#############################
jj                <- which(BP.depth>min.depth & BP.depth<max.depth)

no.IES         <- length(jj)
DDD.list.sub   <- DDD.list[jj]       

for(n in 1:no.IES) {
  
  z.max <- min.depth + max.depth
  nrow  <- nrow(DDD.list.sub[[n]])
  DDD.list.sub[[n]][nrow+1,] <- c(z.max, 0, 0)
  
}



##############
# model SOBP #
##############

resolution         <- 100
plateau.dose.Gy    <- 2

z.max              <- round(z.max, digits=1)
depth.seq          <- seq(from=0, to=z.max, by=z.max/resolution)

min.depth.step     <- ceiling( DDD.list.sub[[1]]$z.cm[which.max(DDD.list.sub[[1]]$D.Gy)] / (z.max/resolution) )
max.depth.step     <- floor( DDD.list.sub[[no.IES]]$z.cm[which.max(DDD.list.sub[[no.IES]]$D.Gy)] / (z.max/resolution) )

primary.weights    <- abs(nlm(f=deviation.proton, p=c(1:no.IES))$estimate)
total.weights      <- round( primary.weights * plateau.dose.Gy / mean(SOBP.dose(primary.weights)[min.depth.step:max.depth.step]), digits=0 )
used.beam.energies <- beam.energies[jj]



############################
# plot SOBP (single field) #
############################

plot1  <- xyplot(SOBP.dose(total.weights) ~ depth.seq,
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