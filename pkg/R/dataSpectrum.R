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
                                   spectrum               = data.frame(particle.no                 = integer(),
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

particles.per.primary <- function(x){
  return(sum(x@spectrum$N.per.primary))
}

Mass.Stopping.Power.MeV.cm2.g <- function(x, stopping.power.source){
  return(  AT.Mass.Stopping.Power( stopping.power.source,
                                   x@spectrum$E.mid.MeV.u,
                                   x@spectrum$particle.no,
                                   AT.material.no.from.material.name(x@target.material))$stopping.power.MeV.cm2.g)
  
}

dose.per.primary <- function(x, stopping.power.source){
  LET.MeV.cm            <- Mass.Stopping.Power.MeV.cm2.g(x, stopping.power.source)
  
  density.g.cm3         <- AT.get.materials.data(AT.material.no.from.material.name(x@target.material))$density.g.cm3
  
  return(sum(LET.MeV.cm * x@spectrum$N.per.primary) / density.g.cm3 * 1.60217657e-10)
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

######################
# Method plot
setMethod(f          = "plot", 
          signature  = c("dataSpectrum"),
          definition = function(x) {
            
            Z <- AT.Z.from.particle.no(x@spectrum$particle.no)$Z
            lattice::xyplot(N.per.primary ~ E.mid.MeV.u,
                            x@spectrum,
                            type     = "S",
                            grid     = TRUE,
                            groups   = Z,
                            ylab     = "particles / Gy",
                            scales   = list(y = list(log = 10)),
                            auto.key = list(title = "Z", 
                                            space = "top", 
                                            columns = max(Z), 
                                            lines = TRUE, 
                                            points = FALSE))
          })