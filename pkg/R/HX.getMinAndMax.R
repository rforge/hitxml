HX.getMinAndMax <- function(beam.spot.grid){
    minMax        <- numeric(4)
	names(minMax) <- c("min.x.mm", "min.y.mm", "max.x.mm", "max.y.mm")
	minMax[1]     <- min(beam.spot.grid$x.mm) - max(beam.spot.grid$focus.X.FWHM.mm) * 2
	minMax[2]     <- min(beam.spot.grid$y.mm) - max(beam.spot.grid$focus.Y.FWHM.mm) * 2
	minMax[3]     <- max(beam.spot.grid$x.mm) + max(beam.spot.grid$focus.X.FWHM.mm) * 2
	minMax[4]     <- max(beam.spot.grid$y.mm) + max(beam.spot.grid$focus.Y.FWHM.mm) * 2
    return(minMax)
}
