HX.compute.field <- function( beam.spot.grid, 
                              field.mm = NULL, 
                              resolution.mm,
                              method = "old"){
    
	# If field size is not given, find reasonable (maximal) size
	if( is.null(field.mm)){
        field.mm    <- HX.getMinAndMax(beam.spot.grid = beam.spot.grid)
	}        

    # Get voxel positions
    x.steps.mm    <- seq(field.mm[1], field.mm[3], by = resolution.mm)
    y.steps.mm    <- seq(field.mm[2], field.mm[4], by = resolution.mm)
    n.x.steps     <- length(x.steps.mm)
    n.y.steps     <- length(y.steps.mm)
    
    # Construct matrix, row 1: x, row 2: y, row 3: fluence
    m             <- matrix(nrow = n.x.steps * n.y.steps, ncol = 3)
    m[,1]         <- rep(x.steps.mm, n.y.steps)
    m[,2]         <- sort(rep(y.steps.mm, n.x.steps))

    m[,3]         <- 0

 
    # Get number of rows to obtain array lengths
    if(method == "old"){
      n.df     <- nrow(beam.spot.grid)
      n.m      <- nrow(m)
      tmp      <- numeric(n.m)
      
     # Feed into C routine to compute fluence as sum of 2D Gauss beam spots
      res <- .C("XML_PBR",
                   n.df = as.integer(n.df),
                   df.x = as.double(beam.spot.grid$x.mm),
                   df.y = as.double(beam.spot.grid$y.mm),
                   df.focusX = as.double(beam.spot.grid$focus.X.FWHM.mm),
                   df.focusY = as.double(beam.spot.grid$focus.Y.FWHM.mm),
                   df.particles = as.double(beam.spot.grid$N.particles),
                   n.m = as.integer(n.m),
                   m.x = as.double(m[,1]),
                   m.y = as.double(m[,2]),
                   m.value = as.double(tmp),
                   PACKAGE = "HITXML")
      
       m[,3] <- res$m.value
     }else{
       # m[,3] <- XML_PBR_new(beam.spot.grid, m)
     }
     return(m)
}