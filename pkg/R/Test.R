if(F){
  
  library(HITXML)
  setwd("D:/svns/svn.R-forge/HITXML/pkg/exec")
  dd <- dataRBE(file.name = "chordom02.rbe")



  ss <- dataSPC(file.name = "12C.H2O.MeV27000.spc")
  sss <- dataSpectrum(ss, 16.5)
    
  dose.per.primary(x = sss, stopping.power.source = "PSTAR")
  
  HX.RBE.LEM(RBE.data = dd,
             Spectrum.data = sss,
             dose.Gy = 1)
}