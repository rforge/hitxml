libamtrack::AT.dose.Gy.from.fluence.cm2(278.29, 6012, 52868500, 1, 3)

df

file.name <- "Test.dat"

  str <- paste(df$energy.MeV.u/1000, df$x.mm, df$y.mm, 
               rowMeans(matrix(data = c(df$focus.X.FWHM.mm, df$focus.Y.FWHM.mm), ncol = 2)),
               df$N.particles)
  write(str,
        file = file.name,
        sep  = "\n")



vv <- 1:9
zz <- 11:19
rowMeans(matrix(data = c(vv, zz), ncol = 2))
