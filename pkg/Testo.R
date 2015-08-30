rm(list = ls())
library(HITXML)
library(lattice)

beam.spot.grid   <- HX.construct.field(field.shape = "square", 
                                       par = 5, 
                                       focus.FWHM.mm = 5.0, 
                                       fluence.cm2 = 100, 
                                       field.par = 100)$beam.spot.grid

a <- system.time(
m <- HX.compute.field(beam.spot.grid, resolution.mm = .5, method = "new"))
b <- system.time(
n <- HX.compute.field(beam.spot.grid, resolution.mm = .5, method = "old"))

a
b


levelplot(m[,3] ~ m[,1]*m[,2],
          useRaster = TRUE,
          aspect    = 1.0,
          cuts      = 100,
          col.regions = grey(0:100/100))
levelplot(n[,3] ~ n[,1]*n[,2],
          useRaster = TRUE,
          aspect    = 1.0,
          cuts      = 100,
          col.regions = grey(0:100/100))

