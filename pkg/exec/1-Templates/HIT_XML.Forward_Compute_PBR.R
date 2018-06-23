############################################
# Script to forward calculate HIT 
# physical beam records 
#
# S. Greilich, Apr. 2013, rev. Sep 2016
# DKFZ Heidelberg, s.greilich@dkfz.de
############################################

rm(list =ls())
# install.packages("HITXML", repos="http://R-Forge.R-project.org", dependencies = "TRUE)
library(HITXML)
library(data.table)
library(ggplot2)


# USER SETTINGS ####
resolution.mm          <- 1
# file.name              <- "HIT-Plan/SOBP.xml"
file.name              <- "ab0001Test.xml"    

# SCRIPT ####
df <- HX.read.PBR(file.name)

if(min(df$focus.X.FWHM.mm) <= 0 || min(df$focus.Y.FWHM.mm) <= 0){
  Stop("PBR contains no focus information (TCU3?), exiting...")
}

n.IES        <- length(unique(df$IES))
  
# Find max field
field.mm    <- HX.getMinAndMax(beam.spot.grid = df) * 1.25
  
  
# Compute all IESs
dt  <-  rbindlist(lapply(1:n.IES,
                         function(i, ddf, f.mm, r.mm){
                           as.data.table( 
                             HX.compute.field( beam.spot.grid = ddf[ddf$IES == i,],
                                               field.mm       = f.mm,
                                               resolution.mm  = r.mm))},
                         ddf  = df,
                         f.mm = field.mm,
                         r.mm = resolution.mm))
names(dt) <- c("x.mm", "y.mm", "fluence.mm2")

# Collapse IESs
dtt  <- dt[,
           .(fluence.cm2 = 100*sum(fluence.mm2)),
           by = .(x.mm, y.mm)]

# Plot field
ii   <- abs(dtt$y.mm) == min(abs(dtt$y.mm))
jj   <- ifelse(dtt$y.mm[which.min(abs(dtt$y.mm))] > 0, 
               dtt$y.mm > 0,
               dtt$y.mm <= 0)
dttt <- dtt[ii&jj,]

ggplot(data = dtt, mapping = aes(x = x.mm, y = y.mm)) +
  theme_bw() +
  geom_tile(aes(fill = fluence.cm2)) +
  geom_hline(aes(yintercept = min(abs(dtt$y.mm))), colour = "white", linetype = 2)


ggplot(data = dttt, mapping = aes(x = x.mm, y = fluence.cm2)) +
  theme_bw() +
  geom_line()

mean(dttt$fluence.cm2[dttt$x.mm > -10 & dttt$x.mm < 10])