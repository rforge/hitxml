read.beam.energy <- function(x, file.names) {
  
  file.name                   <- file.names[x]
  input                       <- scan(file = file.name,
                                      what = "character", strip.white = TRUE, sep = "")
  beam.energy                 <- as.numeric(gsub(",", "", input[grep("energy", input) + 1]))
  
  return(beam.energy)
  
}



get.BP.depth <- function(x, file.names) {

  file.name                   <- file.names[x]
  ddd                         <- read.table( file.name, skip=10 )
  if(ncol(ddd)>2) {
    ddd <- ddd[,-5]
    ddd <- ddd[,-4]
    ddd <- ddd[,-3]
  }
  colnames(ddd)               <- c("z.cm", "dE.dz.MeV.g.cm")
  ddd$D.Gy                    <- (1.602E-10) * ddd$dE.dz.MeV.g.cm
  BP.depth                    <- ddd$z.cm[which.max(ddd$D.Gy)]
  
  return(BP.depth)
  
}



read.DDD <- function(x, file.names) {
  
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



read.SPC <- function(x, file.names) {
  
  file.name              <- file.names[x]
  spc                    <- AT.SPC.read(file.name=file.name, endian="little", flavour="vanilla")$spc
  spc$particle.name      <- AT.particle.name.from.particle.no(particle.no = spc$particle.no)
  spc$Z                  <- AT.Z.from.particle.no(particle.no = spc$particle.no)$Z
  names(spc)[8]          <- paste("dN.dE.MeV.u")
  names(spc)[9]          <- paste("fluence.cm2")
  
  return(spc)
  
}



SOBP.dose <- function(x) {
  
  dose.fun    <- function(idx, z) approx(DDD.list.sub[[idx]]$z.cm, DDD.list.sub[[idx]]$D.Gy, xout=z)$y
      
  for(n in 1:length(jj)) {
    if(n==1) { sum <- x[n] * dose.fun(idx=n, z=depth.seq) }
    else     { sum <- sum + x[n] * dose.fun(idx=n, z=depth.seq) }
  }
  
  return(sum)
  
}



SOBP.LET <- function(x) {
  
  LET.fun    <- function(idx, z) approx(df.total[[idx]]$depth.g.cm2, df.total[[idx]]$LET.keV.um, xout=z)$y
  
  for(n in 1:length(jj)) {
    if(n==1) { sum <- x[n] * LET.fun(idx=n, z=depth.seq) }
    else     { sum <- sum + x[n] * LET.fun(idx=n, z=depth.seq) }
  }
  
  return(sum)
  
}



SOBP.fLET <- function(x) {
  
  fLET.fun    <- function(idx, z) approx(df.total[[idx]]$depth.g.cm2, df.total[[idx]]$fLET.keV.um, xout=z)$y
  
  for(n in 1:length(jj)) {
    if(n==1) { sum <- x[n] * fLET.fun(idx=n, z=depth.seq) }
    else     { sum <- sum + x[n] * fLET.fun(idx=n, z=depth.seq) }
  }
  
  return(sum)
  
}



deviation.carbon <- function(x) {
    
  dev            <- sum( abs(SOBP.dose(abs(x))[min.depth.step:max.depth.step] - mean(SOBP.dose(abs(x))[min.depth.step:max.depth.step])) ) * 1E7
  
  return(dev)
  
}



deviation.proton <- function(x) {
  
  dev            <- sum( abs(SOBP.dose(abs(x))[min.depth.step:max.depth.step] - mean(SOBP.dose(abs(x))[min.depth.step:max.depth.step])) ) * 1E7
  
  return(dev)
  
}