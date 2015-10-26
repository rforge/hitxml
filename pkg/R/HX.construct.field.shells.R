HX.construct.field.shells <- function(d.mm, 
                                      field.size.mm,
                                      focus.FWHM.mm){

	r.seq.mm          <- seq(from = 0, 
	                         to   = field.size.mm/2, 
	                         by   = focus.FWHM.mm)

	get.segments      <- function( radius){
		if(radius == 0.0){return(0)}
	  n.segments <- floor(2 * pi * radius / d.mm)
		if(n.segments <= 0){
		    n.segments <- 1
		}
		phi.seq <- seq( from         = 0, 
					 to           = 2*pi - d.mm / radius, 
					 length.out   = n.segments - 1)
		offset <- 0.5 * (2 * pi - tail(phi.seq, 1))
		return( phi.seq + offset)
	}
	
	phi.seq           <- lapply(r.seq.mm, get.segments)
	df                <- data.frame( r.mm   = rep(r.seq.mm, unlist(lapply(phi.seq, length))),
					                 phi    = unlist(phi.seq))
	
	beam.spot.grid <- data.frame( x.mm            = df$r.mm * sin(df$phi),
								  y.mm            = df$r.mm * cos(df$phi))
	

    return(beam.spot.grid)
}