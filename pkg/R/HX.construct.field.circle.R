HX.construct.field.circle <- function(r.step.mm,
                                      angular.step.mm, 
                                      r.min.mm,
									  field.size.mm,
									  focus.FWHM.mm){

	r.seq.mm          <- seq(from = max(c(r.min.mm - 0.5 * focus.FWHM.mm, 0)), 
	                                to = r.min.mm + field.size.mm + 0.5 * focus.FWHM.mm, 
	                                by = r.step.mm)
	get.segments      <- function( radius, angular.step){
		n.segments <- floor(2 * pi * radius / angular.step)
		if(n.segments <= 0){
		    n.segments <- 1
		}
		phi.seq <- seq( from         = 0, 
					 to           = 2*pi - angular.step / radius, 
					 length.out   = n.segments - 1)
		offset <- 0.5 * (2 * pi - tail(phi.seq, 1))
		return( phi.seq + offset)
	}
	
	phi.seq           <- lapply(r.seq.mm, get.segments, angular.step = angular.step.mm)
	df                <- data.frame( r.mm   = rep(r.seq.mm, unlist(lapply(phi.seq, length))),
					                 phi    = unlist(phi.seq))
	
	beam.spot.grid <- data.frame( x.mm            = df$r.mm * sin(df$phi),
								  y.mm            = df$r.mm * cos(df$phi))
	

    return(beam.spot.grid)
}