fluence.cm2.mid.SOBP <- total[which.min(abs(LET.depths.g.cm2 - 20.0))]
fluence.cm2.mid.SOBP.requested <- 2.8e6

scale.factor <- fluence.cm2.mid.SOBP.requested / fluence.cm2.mid.SOBP