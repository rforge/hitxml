writeDDD <- function(ddd, projectile, energy.MeV){
  output <- c("!filetype    ddd", 
              "!fileversion    19980520", 
              paste0("!filedate    ", date()), 
              paste0("!projectile    ", projectile),
              "!material      H2O",
              "!composition   H2O",
              "!density 1",
              paste0("!energy ", energy.MeV),
              "   z[g/cm**2] dE/dz[MeV/(g/cm**2)] FWHM1[g/cm**2] factor FWHM2[g/cm**2]\n!ddd\n")
  write(output, "test.ddd", sep = "\n")
}


writeDDD(NULL, "1H", "23.66")