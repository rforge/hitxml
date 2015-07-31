empty.spectrum <- function(nrow){
  return(matrix(nrow     = nrow,
                ncol     = 4,
                dimnames = list(NULL,
                                c("particle.no",
                                  "E.MeV.u",
                                  "dE.MeV.u",
                                  "N"))))
}
################################
# dataSpectrum CLASS
################################
setClass( Class            = "dataSpectrum",
          slots            = c( depth.g.cm2            = "numeric",
                                spectrum               = "matrix"),
          prototype        = list( depth.g.cm2            = numeric(),
                                   spectrum               = empty.spectrum(0)) )
################################
# Constructor
################################
dataSpectrum <- function(SPC.data, depth.g.cm2){
  res   <- lapply(1:length(depth.g.cm2),
                  function(i, s, d){
                    new("dataSpectrum",
                        depth.g.cm2         = d[i],
                        spectrum            = spectrum.at.depth.g.cm2(s@spectra, d[i]))},
                  s = SPC.data,
                  d = depth.g.cm2)
    
}

particles <- function(x){
  return(sum(x@spectrum[,"N"]))
}

Mass.Stopping.Power.MeV.cm2.g <- function(x, stopping.power.source, target.material){
  return(  AT.Mass.Stopping.Power( stopping.power.source,
                                   x@spectrum[,"E.MeV.u"],
                                   x@spectrum[,"particle.no"],
                                   AT.material.no.from.material.name(target.material))$stopping.power.MeV.cm2.g)
  
}

dose.from.spectrum.Gy <- function(x, stopping.power.source, target.material){
  LET.MeV.cm            <- Mass.Stopping.Power.MeV.cm2.g(x, stopping.power.source)
  
  density.g.cm3         <- AT.get.materials.data(AT.material.no.from.material.name(target.material))$density.g.cm3
  
  return(sum(LET.MeV.cm * x@spectrum$N) / density.g.cm3 * 1.60217657e-10)
}

###########################################
# R function for interpolation of spc files
# with depth, was in libamtrack
# moved to HITXML Jul15, sgre
###########################################

spectrum.at.depth.g.cm2 <- function(spectra, depth.g.cm2, interpolate = TRUE)
{
  depth.step.of.spc      <- unique(spectra[,"depth.step"])
  depth.g.cm2.of.spc     <- unique(spectra[,"depth.g.cm2"])
  depth.step.interp      <- approx( x    = depth.g.cm2.of.spc,
                                    y    = depth.step.of.spc,
                                    xout = depth.g.cm2,
                                    rule = 2)$y
  
  depth.step.int         <- floor(depth.step.interp)
  depth.step.frac        <- depth.step.interp - depth.step.int
  spectrum.upstream      <- spectra[spectra[,"depth.step"] == depth.step.int,] 
  spectrum.downstream    <- spectra[spectra[,"depth.step"] == (depth.step.int+1),]  
  spectrum.interp                 <- empty.spectrum(nrow(spectrum.upstream))
  spectrum.interp[,"E.MeV.u"]     <- spectrum.upstream[,"E.MeV.u"]
  spectrum.interp[,"dE.MeV.u"]    <- spectrum.upstream[,"dE.MeV.u"]
  spectrum.interp[,"particle.no"] <- spectrum.upstream[,"particle.no"]
  if(interpolate){
    spectrum.interp[,"N"]           <- (1 - depth.step.frac) * spectrum.upstream[,"N.per.primary"] + depth.step.frac * spectrum.downstream[,"N.per.primary"]
  }else{
    spectrum.interp[,"N"]           <- spectrum.upstream[,"N.per.primary"]
  }
  return(spectrum.interp)
}

######################
# Method plot
setMethod(f          = "plot", 
          signature  = c("dataSpectrum"),
          definition = function(x) {
            
            Z <- AT.Z.from.particle.no(x@spectrum[,"particle.no"])$Z
            lattice::xyplot(x@spectrum[,"N"] ~ x@spectrum[,"E.MeV.u"],
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



setMethod(f          = "+", 
          signature  = c("dataSpectrum", "dataSpectrum"),
          definition = function(e1, e2) {
            
          new.spectrum                <- e1@spectrum
          new.spectrum[,"N"]          <- e1@spectrum[,"N"] + e2@spectrum[,"N"]
              
          return(new("dataSpectrum",
                     depth.g.cm2         = e1@depth.g.cm2,
                     spectrum            = new.spectrum))
})

setMethod(f          = "*", 
          signature  = c("dataSpectrum", "numeric"),
          definition = function(e1, e2) {
            
            new.spectrum                <- e1@spectrum
            new.spectrum[,"N"]              <- new.spectrum[,"N"] * e2

            new("dataSpectrum",
                depth.g.cm2         = e1@depth.g.cm2,
                spectrum            = new.spectrum)
          })
