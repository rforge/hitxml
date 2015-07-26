################################
# dataSpectrum CLASS
################################
setClass( Class            = "dataSpectrum",
          slots            = c( projectile             = "character",
                                beam.energy.MeV.u      = "numeric",
                                target.material        = "character",
                                peak.position.g.cm2    = "numeric",
                                depth.g.cm2            = "numeric",
                                spectrum               = "data.frame"),
          prototype        = list( projectile             = character(),
                                   beam.energy.MeV.u      = numeric(),
                                   target.material        = character(),
                                   peak.position.g.cm2    = numeric(),
                                   depth.g.cm2            = numeric(),
                                   spectra                = data.frame(particle.no                 = integer(),
                                                                       E.low.MeV.u                 = numeric(),
                                                                       E.mid.MeV.u                 = numeric(),
                                                                       E.high.MeV.u                = numeric(),
                                                                       dE.MeV.u                    = numeric(),
                                                                       N.per.primary               = numeric())) )
################################
# Constructor
################################
dataSpectrum <- function(SPC.data, depth.g.cm2){
  res   <- SPC.spectrum.at.depth.g.cm2(SPC.data@spectra, depth.g.cm2)
    
  new("dataSpectrum",
      projectile          = SPC.data@projectile,
      beam.energy.MeV.u   = SPC.data@beam.energy.MeV.u,
      target.material     = SPC.data@target.material,
      peak.position.g.cm2 = SPC.data@peak.position.g.cm2,
      depth.g.cm2         = depth.g.cm2,
      spectrum            = res)
}

###########################################
# R function for interpolation of spc files
# with depth, was in libamtrack
# moved to HITXML Jul15, sgre
###########################################
SPC.spectrum.at.depth.g.cm2 <- function(spc, depth.g.cm2, interpolate = TRUE)
{
  depth.step.of.spc      <- unique(spc$depth.step)
  depth.g.cm2.of.spc     <- unique(spc$depth.g.cm2)
  depth.step.interp      <- approx( x    = depth.g.cm2.of.spc,
                                    y    = depth.step.of.spc,
                                    xout = depth.g.cm2,
                                    rule = 2)$y
  
  depth.step.int         <- floor(depth.step.interp)
  depth.step.frac        <- depth.step.interp - depth.step.int
  spc.before             <- spc[spc$depth.step == depth.step.int,] 
  if(interpolate){
    spc.after              <- spc[spc$depth.step == (depth.step.int+1),]  
    
    spc.interp             <- spc.before
    spc.interp$depth.step  <- depth.step.interp
    spc.interp$depth.g.cm2 <- depth.g.cm2
    spc.interp$N.per.primary <- (1 - depth.step.frac) * spc.before$N.per.primary + depth.step.frac * spc.after$N.per.primary
    
    return(spc.interp)
  }else{
    return(spc.before)
  }
}
