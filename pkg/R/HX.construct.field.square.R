HX.construct.field.square <- function(d.mm, 
                                      field.size.mm,
									                    focus.FWHM.mm){

	## Construct beam spot position grid
	# The minimum side length is homogenous field + focus width
  # Always use odd number of spots (to create overlap in beam's eye view)
	min.grid.size.mm   <- field.size.mm + 2 * focus.FWHM.mm
  n.steps            <- floor(min.grid.size.mm / d.mm)
  if(n.steps%%2 == 1){
    n.steps <- n.steps - 1
  } 
  grid.size.mm       <- n.steps * d.mm
  steps.mm           <- seq(from         = -grid.size.mm/2,
                            to           = grid.size.mm/2,
                            by           = d.mm) 
    
    # Dataframe holds beam spot positions, focus size and particle numbers
    # It anticipates the structure of the HIT XML plan / record files
    beam.spot.grid  <- expand.grid( x.mm            = steps.mm,
                                    y.mm            = steps.mm)

    # construct aux vectors to sort beam spots to minimize jumps, x.scan.direction for illustration only
    x.scan.direction <- rep( c(rep("+", length(steps.mm)), rep("-", length(steps.mm))), 
							length.out = nrow(beam.spot.grid))
    spot.line.no     <- rep( c(1:length(steps.mm), length(steps.mm):1),
							 length.out = nrow(beam.spot.grid))
    spot.no          <- sort(rep( 1:length(steps.mm), length(steps.mm)))
 
 	beam.spot.grid   <- beam.spot.grid[order(spot.no, spot.line.no),]
    beam.spot.grid$spot.no <- 1:nrow(beam.spot.grid)
 

    return(beam.spot.grid)
}