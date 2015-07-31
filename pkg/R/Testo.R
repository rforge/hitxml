if(F){
  rm(list = ls())

library(HITXML)

spc.path <- "D:/04 - Risoe, DKFZ/03 - Methodik/11-20/20 - TRiP/04 - TRiP Basic Data/HIT/03 - TRiP98DATA_HIT-20131120/SPC/12C/RF3MM"

spc <- dataSPC(file.name = file.path(spc.path, "FLUKA_NEW3_12C.H2O.MeV27000.spc"), 
               endian = "little",
               redistribute = FALSE)
sspc <- redistribute.spc(spc)

s   <- dataSpectrum(spc, 10.0)
ss  <- dataSpectrum(sspc, 10.0)

plot(s[[1]])
plot(ss[[1]])


}