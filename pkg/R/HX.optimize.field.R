HX.optimize.field <- function( par.start,
                               focus.FWHM.mm,
							   fluence.cm2,
                               N.min,
							   field.shape,
							   field.par,
							   n.IES,
							   resolution.mm = 1,
							   spot.distance.mm,
							   plot = FALSE){

	# Do actual optimization using cost.function
	# DEBUG: d.mm <- d.start.mm
  # If positive value for spot distance is given, don't optimize
	if(spot.distance.mm <= 0){
	  res <- optim( par            = par.start,
	              fn             = HX.cost.function,
	              gr             = NULL,
	              field.shape    = field.shape,
				        focus.FWHM.mm  = focus.FWHM.mm,
	              fluence.cm2    = fluence.cm2,
	              field.par      = field.par,
	              N.min          = N.min,
				        resolution.mm  = resolution.mm,
			          method         = "L-BFGS-B",
	              lower          = 1,
	              upper          = 2 * focus.FWHM.mm,
			          control        = list( factr = 1e9))
  	par             <- res$par
	}else{
	  cat("Using fixed spot distance of ", spot.distance.mm, " mm\n")
	  par             <- c(spot.distance.mm, 0) 
	}
  
	optim.field     <- HX.construct.field(field.shape                = field.shape,
	                                      par                        = par, 
	                                      focus.FWHM.mm              = focus.FWHM.mm, 
	                                      fluence.cm2                = fluence.cm2, 
	                                      field.par                  = field.par,
									      compute.with.resolution.mm = resolution.mm)
	switch( match(field.shape, c("square", "circle")),
				cat("Optimized field with spot distance d = ", par[1], "mm.\n"),
				cat("Optimized field with radial spot distance = ", par[1], "mm and angular distance =", par[2], "mm.\n"),
				cat("wrong field shape!"))
	
 	# TODO: apply travelling salesman problem to find best scan path
#	new.order      <- as.integer( solve_TSP( TSP( dist( optim.field$beam.spot.grid[c("x.mm", "y.mm")], 
#	                                                    method = c("euclidean", "manhattan")[2]))))
#	optim.field$beam.spot.grid <- optim.field$beam.spot.grid[new.order,]
#	optim.field$beam.spot.grid$spot.no <- 1:nrow(optim.field$beam.spot.grid)
 
	if(plot == TRUE){
        require(lattice)
        
 		m               <- HX.compute.field( beam.spot.grid = optim.field$beam.spot.grid, 
 	    	                                 field.mm       = NULL,
											 resolution.mm  = resolution.mm)
        # Plot spot positions
        custom.panel    <- function(spot.no = spot.no, x = x, y = y, ...){
        	panel.xyplot(x = x, y = y, ...)
		    for(i in 1:length(spot.no)){
        		panel.text( spot.no[i],
					    x = x[i],
					    y = y[i],
					    cex = 10 / sqrt(nrow(optim.field$beam.spot.grid)))
        	}
        }
        
        aspect <- (max(optim.field$beam.spot.grid$y.mm) - min(optim.field$beam.spot.grid$y.mm)) /
                     (max(optim.field$beam.spot.grid$x.mm) - min(optim.field$beam.spot.grid$x.mm))

        plot(
        	xyplot( y.mm ~ x.mm,
				optim.field$beam.spot.grid,
				type = 'p',
				pch  = 3,
				cex  = 10 / sqrt(nrow(optim.field$beam.spot.grid)),
				spot.no = optim.field$beam.spot.grid$spot.no,
				panel = custom.panel,
				xlim = c(min(m[,1]), max(m[,1])),
				ylim = c(min(m[,2]), max(m[,2])),
			    main = paste('Optimized plan - entire field (IES no. ', n.IES, ')', sep = ""), 
			    sub  = 'beam spot positions and numbers',
				aspect = abs(aspect))
		)		
		# Plot full field
         
		y.line.mm       <- min(abs(m[,2]))
		custom.panel    <- function(...){
            panel.levelplot(...)
			panel.xyplot( optim.field$beam.spot.grid$x.mm, 
				          optim.field$beam.spot.grid$y.mm,
				          type = 'p',
				 		  pch  = 3)
		    panel.abline(h = y.line.mm, lty = 2)
		}
		
		plot(
		levelplot(m[,3]*100 ~ m[,1]*m[,2],
			      panel     = custom.panel,
				  cex       = 20 / sqrt(nrow(optim.field$beam.spot.grid)),
			      xlab      = 'x / mm',
			      ylab      = 'y / mm',
			      main      = paste('Optimized plan - entire field (IES no. ', n.IES, ')', sep = ""),
			      sub       = 'fluence / (1/cm2)',
			      aspect    = aspect)
        )
		
		# Plot cross section
 		plot(
		xyplot(m[,3]*100 ~ m[,1],
		       type   = 's',
 		       subset = abs(m[,2]) == y.line.mm,
		       xlab   = 'x / mm',
		       ylab   = 'fluence / (1/cm2)',
		       main   = paste('Optimized plan - cross section (IES no. ', n.IES, ')', sep = ""))
        )
	    
		# Plot homogenous part, with relative deviation from requested fluence
 		if(field.shape == "square" & length(field.par) == 1){
 		  field.par <- c(field.par, field.par)
 		}
 		m               <- HX.compute.field( beam.spot.grid = optim.field$beam.spot.grid, 
 	    	                                 field.mm       = c( -field.par[1]/2, -field.par[2]/2, field.par[1]/2, field.par[2]/2),
											 resolution.mm  = resolution.mm)
 		
 		if(fluence.cm2 == 0){
 		  stop("Problem with zero fluence!")
 		}
 		
 		plot(
		levelplot( (m[,3]*100 - fluence.cm2) /fluence.cm2 * 100 ~ m[,1]*m[,2],
			      panel = custom.panel,
				  xlab  = 'x / mm',
			      ylab  = 'y / mm',
			      main      = paste('Optimized plan - homogenous central area (IES no. ', n.IES, ')', sep = ""),
			      sub       = 'deviation from requested fluence / %'),
		)
 		
	}
    return(optim.field)
}
