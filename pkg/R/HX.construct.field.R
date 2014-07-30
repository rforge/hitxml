HX.construct.field <- function(field.shape,
                               par, 
                               focus.FWHM.mm, 
                               fluence.cm2, 
                               field.par,
							   compute.with.resolution.mm = NULL){
    
    switch( match(field.shape, c("square", "circle")),
	        beam.spot.grid <- HX.construct.field.square( d.mm                       = par[1], 
		                                                 field.size.mm              = field.par[1],
									                     focus.FWHM.mm              = focus.FWHM.mm),
	        beam.spot.grid <- HX.construct.field.circle( r.step.mm                  = par[1], 
		                                                 angular.step.mm            = par[2],
	        											 field.size.mm              = field.par[1],
	        											 r.min.mm                   = field.par[2],
									                     focus.FWHM.mm              = focus.FWHM.mm),
            cat("wrong field shape!"))
    
    beam.spot.grid$focus.X.FWHM.mm <- focus.FWHM.mm
    beam.spot.grid$focus.Y.FWHM.mm <- focus.FWHM.mm
    beam.spot.grid$N.particles     <- 1
    
    
    mean.fluence.mm2              <- 0
	sd.fluence.mm2                <- 0
	N.particles.total             <- 0

    # If given, compute field with given resolution and scale total fluence
    if( !is.null(compute.with.resolution.mm)){
    	m                  <- HX.compute.field(beam.spot.grid = beam.spot.grid, 
                                               resolution.mm  = compute.with.resolution.mm)
		switch( match(field.shape, c("square", "circle")),
				ii <- (m[,1] >= -field.par[1]/2) &
				      (m[,1] <=  field.par[1]/2) &
					  (m[,2] >= -field.par[1]/2) &
					  (m[,2] <=  field.par[1]/2),
				ii <- (sqrt(m[,1]^2+m[,2]^2) <= field.par[1] + field.par[2]) &
				      (sqrt(m[,1]^2+m[,2]^2) >= field.par[2]),
				cat("wrong field shape!"))
        mean.fluence.mm2              <- mean(m[ii,3])
        sd.fluence.mm2                <- sd(m[ii,3]) / mean.fluence.mm2
        N.particles.total             <- floor((fluence.cm2/100) / mean.fluence.mm2) + 1
        beam.spot.grid$N.particles    <- N.particles.total * beam.spot.grid$N.particles
    }

    return(list( beam.spot.grid    = beam.spot.grid,
                 N.particles.total = N.particles.total,
                 mean.fluence.mm2  = mean.fluence.mm2 * N.particles.total,
                 sd.fluence.mm2    = sd.fluence.mm2))
}