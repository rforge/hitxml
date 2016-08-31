library(lattice)

# Run HIT_XML.Create_SOBP.R
depths.cm <- sapply(eff.spectra.at.depth, function(x){x@depth.g.cm2})

which(depths.cm == 1.25)
which(depths.cm == 7.75)

entr.spc <- as.data.frame(eff.spectra.at.depth[[6]]@spectrum)
sobp.spc <- as.data.frame(eff.spectra.at.depth[[32]]@spectrum)

entr.spc$location <- "entrance (1 cm)"
sobp.spc$location <- "sobp (7.5 cm)"

df    <- rbind(entr.spc, sobp.spc)
df$particle.name <- factor(df$particle.no, 
                           levels = c(1002,2004,3006,4008,5010,6012),
                           labels = c("H", "He", "Li", "Be", "B", "C"))
df$LET.keV.um    <- AT.Stopping.Power("ICRU", df$E.MeV.u, df$particle.no, 1)$stopping.power.keV.um

xyplot(N/10 ~ E.MeV.u|location,
       df,
       type   = "l",
       groups = particle.name,
       grid   = TRUE,
       as.table = TRUE,
       auto.key  = list(space = "top", columns = 6, points = FALSE, lines = TRUE),
       scales = list(log = TRUE, equispaced.log = FALSE),
       ylim   = c(9e1, 4e5),
       xlab   = "energy / (MeV/u)",
       ylab   = "fluence / cm-2")


ii <- df$location == "entrance (1 cm)"
hh <- get.multivariate.weighted.histogram( df$LET.keV.um[ii],
                                           rep(1, sum(ii)),
                                           df$N[ii],
                                           0.1,
                                           200,
                                           100,
                                           TRUE)
gg <- get.multivariate.weighted.histogram( df$LET.keV.um[!ii],
                                           rep(1, sum(!ii)),
                                           df$N[!ii],
                                           0.1,
                                           200,
                                           100,
                                           TRUE)
hh$location <- "entrance (1 cm)"
gg$location <- "sobp (7.5 cm)"
dff <- rbind(hh,gg)

xyplot(frequency ~ x|location,
       dff,
       type   = "s",
       grid   = TRUE,
       as.table = TRUE,
       scales = list(log = TRUE, equispaced.log = FALSE),
       ylim   = c(8e3, 5e7),
       xlab   = "LET / (keV/um)",
       ylab   = "frequency / a.u.")
