HX.read.DDD <- function(x, file.names) {
  
  file.name <- DDD.file.names[1]
  
  file.name                   <- file.names[x]
  ddd                         <- read.table( file.name, skip=10 )
  if(ncol(ddd)>2) {
    ddd <- ddd[,-5]
    ddd <- ddd[,-4]
    ddd <- ddd[,-3]
  }
  colnames(ddd)               <- c("z.cm", "dE.dz.MeV.g.cm")
  ddd$D.Gy                    <- (1.602E-10) * ddd$dE.dz.MeV.g.cm
  
  return(ddd)
  
}