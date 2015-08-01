#' @title Create empty spectrum
#' 
#' @description \code{spectrum} return an
#' empty spectrum, i.e. a matrix with four
#' columns particle, kinetic energy (mid bin, in MeV.u),
#' bin width and number of particles in the mid
#' (absolute, not density)
#' 
#' @param nrow Number of rows in the spectrum matrix
#' 
#' @seealso \code{\link{dataSpectrum}}
#' 
#' @family Data structures
#' 
#' @author Steffen Greilich
spectrum <- function(nrow){
  return(matrix(nrow     = nrow,
                ncol     = 4,
                dimnames = list(NULL,
                                c("particle.no",
                                  "E.MeV.u",
                                  "dE.MeV.u",
                                  "N"))))
}

#' @title Class for particle spectra
#' 
#' @description \code{dataSpectrum} is the main class
#' to hold particle spectra information, i.e. the number of
#' particles differential in charge and energy
#' 
#' @slot depth.g.cm2 Depth at which the spectrum is located in g/cm2
#' @slot spectrum \code{\link{spectrum}} (i.e. matrix) structure containing the actual particle numbers 
#'
#' @family Data structures
#' 
#' @author Steffen Greilich

setClass( Class            = "dataSpectrum",
          slots            = c( depth.g.cm2            = "numeric",
                                spectrum               = "matrix"),
          prototype        = list( depth.g.cm2            = numeric(),
                                   spectrum               = spectrum(0)) )

#' @title Constructor for \code{\link{dataSpectrum}} class
#' 
#' @description Constructs a dataSpectrum object from a given dataSPC object, 
#' returning the spectrum at a specific depth
#' 
#' @details Interpolation is done by linear mixture of the spectra which can lead
#' to funny shapes such as double peaks. This is, however, the only way without
#' using explicit transport calculation and also employed for example in TRiP98.
#' 
#' @param SPC.data dataSPC object holding spectrum information for multiple sample depths
#' @param depth.g.cm2 Depth at which the spectrum is located in g/cm2. This can be any
#' depth, spectra will be interpolated between the adjacent sample depths.
#' 
#' @family Data structures
#' 
#' @author Steffen Greilich
dataSpectrum <- function(SPC.data, depth.g.cm2){
  res   <- lapply(1:length(depth.g.cm2),
                  function(i, s, d){
                    new("dataSpectrum",
                        depth.g.cm2         = d[i],
                        spectrum            = spectrum.at.depth.g.cm2(s@spectra, d[i]))},
                  s = SPC.data,
                  d = depth.g.cm2)
    
}

#' @title Total number of particles in a spectrum
#' 
#' @description Return sum of particles
#' 
#' @param x Object of class \code{\link{dataSpectrum}} or list of objects of this class
total.n.particles <- function(x){
  FUN <- function(xx){sum(xx@spectrum[,"N"])}
  if(class(x) == "dataSpectrum"){
    return(FUN(x))
  }else{
    return(sapply(x, FUN))
  }
}

#' @title Mass stopping power for all energy bins in a spectrum
#' 
#' @description (List of) vector(s) of mass stopping power values (in the order of particle number and energies)
#' of the given spectrum/a.
#' 
#' @details Uses \code{\link[libamtrack]{libamtrack}} stopping power function
#' 
#' @param x Object of class \code{\link{dataSpectrum}} or list of objects of this class
#' @param stopping.power.source Descriptor for source of stopping power data (\code{\link[libamtrack]{stopping.power.source}})
#' @param target.material Descriptor for target material (\code{\link[libamtrack]{material.no}})
Mass.Stopping.Power.MeV.cm2.g <- function(x, stopping.power.source, target.material){
  FUN <- function(xx, s, t){ AT.Mass.Stopping.Power( s,
                                                     xx@spectrum[,"E.MeV.u"],
                                                     xx@spectrum[,"particle.no"],
                                                     AT.material.no.from.material.name(t))$stopping.power.MeV.cm2.g}
  
  if(class(x) == "dataSpectrum"){
    return( FUN(x,
                s = stopping.power.source,
                t = target.material))
  }else{ return(  lapply(x,
                  FUN,
                  s = stopping.power.source,
                  t = target.material))
  }
  
}

#' @title Dose for a spectrum
#' 
#' @description (Vector of) dose in Gy
#' 
#' @details Uses \code{\link[libamtrack]{libamtrack}} stopping power function
#' 
#' @param x Object of class \code{\link{dataSpectrum}} or list of objects of this class
#' @param stopping.power.source Descriptor for source of stopping power data (\code{\link[libamtrack]{stopping.power.source}})
#' @param target.material Descriptor for target material (\code{\link[libamtrack]{material.no}})
dose.Gy <- function(x, stopping.power.source, target.material){
  FUN <- function(xx, s, t, d){ m <- Mass.Stopping.Power.MeV.cm2.g(xx, s, t)
                                sum(m * xx@spectrum[,"N"]) / d * 1.60217657e-10}
  
  density.g.cm3 <- AT.get.materials.data(AT.material.no.from.material.name(target.material))$density.g.cm3

  if(class(x) == "dataSpectrum"){
    return(FUN(x,
               s = stopping.power.source,
               t = target.material,
               d = density.g.cm3))
  }else{
    return(sapply(x,
                  FUN,
                  s = stopping.power.source,
                  t = target.material,
                  d = density.g.cm3))
  }
}

#' R function for interpolation of spc files
#' with depth, was in libamtrack and moved to HITXML Jul15, sgre
spectrum.at.depth.g.cm2 <- function(spectra, depth.g.cm2, interpolate = TRUE)
{
  depth.step.of.spc      <- sort(unique(spectra[,"depth.step"]))
  depth.g.cm2.of.spc     <- sort(unique(spectra[,"depth.g.cm2"]))
  depth.step.interp      <- approx( x    = depth.g.cm2.of.spc,
                                    y    = depth.step.of.spc,
                                    xout = depth.g.cm2,
                                    rule = 2)$y
  
  depth.step.int         <- floor(depth.step.interp)
  depth.step.frac        <- depth.step.interp - depth.step.int
  spectrum.upstream      <- spectra[spectra[,"depth.step"] == depth.step.int,] 
  spectrum.interp                 <- spectrum(nrow(spectrum.upstream))
  spectrum.interp[,"E.MeV.u"]     <- spectrum.upstream[,"E.MeV.u"]
  spectrum.interp[,"dE.MeV.u"]    <- spectrum.upstream[,"dE.MeV.u"]
  spectrum.interp[,"particle.no"] <- spectrum.upstream[,"particle.no"]
  if(interpolate & depth.step.frac != 0){
    spectrum.downstream    <- spectra[spectra[,"depth.step"] == (depth.step.int+1),]  
    spectrum.interp[,"N"]           <- (1 - depth.step.frac) * spectrum.upstream[,"N.per.primary"] + depth.step.frac * spectrum.downstream[,"N.per.primary"]
  }else{
    spectrum.interp[,"N"]           <- spectrum.upstream[,"N.per.primary"]
  }
  return(spectrum.interp)
}

#' Method plot
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
