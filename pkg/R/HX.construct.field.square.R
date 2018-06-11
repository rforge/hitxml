HX.construct.field.square <- function(d.mm, 
                                      field.size.mm,
									                    focus.FWHM.mm){

	## Construct beam spot position grid
	# The minimum side length is homogenous field + focus width
  # Always use odd number of spots (to create overlap in beam's eye view)
	min.grid.size.mm   <- field.size.mm + 2 * focus.FWHM.mm
  n.steps            <- floor(min.grid.size.mm / d.mm)
  ii                 <- n.steps%%2 == 1
  n.steps[ii]        <- n.steps[ii] - 1
  grid.size.mm       <- n.steps * d.mm

  steps.x.mm        <- seq(from         = -grid.size.mm[1]/2,
                           to           = grid.size.mm[1]/2,
                           by           = d.mm) 
  if(length(grid.size.mm) == 1){
    steps.y.mm <- steps.x.mm
  }else{
    steps.y.mm        <- seq(from         = -grid.size.mm[2]/2,
                             to           = grid.size.mm[2]/2,
                             by           = d.mm) 
  }
  # Dataframe holds beam spot positions, focus size and particle numbers
  # It anticipates the structure of the HIT XML plan / record files
  beam.spot.grid  <- expand.grid( x.mm            = steps.x.mm,
                                  y.mm            = steps.y.mm)

  # construct aux vectors to sort beam spots to minimize jumps, x.scan.direction for illustration only
  x.scan.direction <- rep( c(rep("+", length(steps.x.mm)), rep("-", length(steps.x.mm))), 
						length.out = nrow(beam.spot.grid))
  spot.line.no     <- rep( c(1:length(steps.x.mm), length(steps.x.mm):1),
						 length.out = nrow(beam.spot.grid))
  spot.no          <- sort(rep( 1:length(steps.y.mm), length(steps.x.mm)))
 
  beam.spot.grid   <- beam.spot.grid[order(spot.no, spot.line.no),]
  beam.spot.grid$spot.no <- 1:nrow(beam.spot.grid)
 

  return(beam.spot.grid)
}