HX.cost.function <- function(field.shape,
                             par,
							               focus.FWHM.mm, 
                             fluence.cm2, 
                             field.par,
							               N.min,
							               resolution.mm){

	# Construct field with current spot distance
	cur.field <- HX.construct.field(   field.shape                = field.shape,
                                     par                        = par, 
                                     focus.FWHM.mm              = focus.FWHM.mm, 
                                     fluence.cm2                = fluence.cm2, 
                                     field.par                  = field.par,
							                       compute.with.resolution.mm = resolution.mm)
    
	# Compute fluence variance
    cost       <- cur.field$sd.fluence.mm2^2 
    
    # Penalize if less particles then minimum machine limit
 	if(cur.field$N.particles.total < N.min){
        cost <- cost + (N.min - cur.field$N.particles.total)^4
    }

    # Output
	par.out <- ""
	for(i in 1:length(par)){
	   par.out <- paste( par.out, 
	                     par[i],
						 "| ",
						 sep = "")
	}
	cat( par.out,
	     "N", sprintf("%4.3e", cur.field$N.particles.total),
         "| sd", sprintf("%4.3e", cur.field$sd.fluence.mm2),
         "| cost", sprintf("%4.3e", cost),
	     "\n")

    return(cost)
}
