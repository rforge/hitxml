if(F){
  
  rm(list = ls())
  library(HITXML)
  setwd("D:/svns/svn.R-forge/HITXML/pkg/exec")
  sss <- dataSPCset(spc.path = "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM")
  ss <- get.spc(SPC.set = sss, beam.energy.MeV.u = 105)


  }