if(F){
  
  dd <- dataRBE(file.name = "chordom02.rbe")

beta.X(dd)

particle.no <- c("1H", "1H", "12C", "14U")
E.MeV.u     <- c(1.0, 10.0, 200.0, 200.0)

RBE.initial(x = dd, projectile = particle.no, E.MeV.u)
alpha.ion(x = dd, projectile = particle.no, E.MeV.u)
beta.ion(x = dd, projectile = particle.no, E.MeV.u)

  ss <- dataSPC(file.name = "12C.H2O.MeV27000.spc")
  sss <- dataSpectrum(ss, 10.5)
}