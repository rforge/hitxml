if(F){
  
  library(HITXML)
  setwd("D:/svns/svn.R-forge/HITXML/pkg/exec")
  dd <- dataRBE(file.name = "chordom02.rbe")


  ss <- dataSPC(file.name = "12C.H2O.MeV27000.spc")
  DDD <- depth.dose(ss)

  DDD$D.Gy <- DDD$D.Gy / DDD$D.Gy[1]
  xyplot(D.Gy ~ depth.g.cm2,
         DDD)

  for(i in 1:nrow(DDD)){
    tmp <- HX.RBE.LEM(RBE.data = dd,
                     Spectrum.data = dataSpectrum(ss, DDD$depth.g.cm2[i]),
                     dose.Gy = DDD$D.Gy[i])
    DDD$RBE[i] <- tmp[1]
    DDD$D.biol.Gy[i] <- tmp[3]
    cat("Did ", i, "\n")
  }
  
  
}