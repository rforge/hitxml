HX.read.beam.energy <- function(x, file.names) {
  
  file.name                   <- file.names[x]
  input                       <- scan(file = file.name,
                                      what = "character", strip.white = TRUE, sep = "")
  beam.energy                 <- as.numeric(gsub(",", "", input[grep("energy", input) + 1]))
  
  return(beam.energy)
  
}
