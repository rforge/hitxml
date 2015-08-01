if(F){
  
  rm(list = ls())
  library(HITXML)
  setwd("D:/svns/svn.R-forge/HITXML/pkg/exec")
  sss <- dataSPCset(spc.path = "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM")
  ss <- get.spc(SPC.set = sss, beam.energy.MeV.u = 105)
  s  <- dataSpectrum(ss, depth.g.cm2 = 2.0)
  plot(s)

  
  
  rm(list = ls())
  library(HITXML)
  setwd("D:/svns/svn.R-forge/HITXML/pkg/exec")
  ddd <- dataDDDset(ddd.path = "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/DDD/12C/RF3MM")
  dd <- get.ddd(DDD.set = ddd, beam.energy.MeV.u = 105)
  plot(dd)
  
  
    }