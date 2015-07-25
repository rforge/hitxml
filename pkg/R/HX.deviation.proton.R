HX.deviation.proton <- function(x) {
  
  dev            <- sum( abs(HX.SOBP.dose(abs(x))[min.depth.step:max.depth.step] - mean(HX.SOBP.dose(abs(x))[min.depth.step:max.depth.step])) ) * 1E7
  
  return(dev)
  
}