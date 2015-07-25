HX.read.SPC <- function(x, file.names) {
  
  file.name              <- file.names[x]
  spc                    <- AT.SPC.read(file.name=file.name, endian="little", flavour="vanilla")$spc
  spc$particle.name      <- AT.particle.name.from.particle.no(particle.no = spc$particle.no)
  spc$Z                  <- AT.Z.from.particle.no(particle.no = spc$particle.no)$Z
  names(spc)[8]          <- paste("dN.dE.MeV.u")
  names(spc)[9]          <- paste("fluence.cm2")
  
  return(spc)
  
}
