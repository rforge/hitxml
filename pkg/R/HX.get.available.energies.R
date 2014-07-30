HX.get.available.energies <- function(energy.chosen.MeV.u, df.libc){
  
  # if user value is not in libC suggest closest values
  energy.diff        <- abs(energy.chosen.MeV.u - df.libc$energy.MeV.u)
  closest.energy.idx <- which(min(energy.diff) == energy.diff)
  neighbor.idx       <- -1:1
  if(closest.energy.idx == 1){                    neighbor.idx <-  0:1}
  if(closest.energy.idx == length(energy.diff)){  neighbor.idx <- -1:0}
  return(df.libc$energy.MeV.u[closest.energy.idx + neighbor.idx])
}