HX.deviation.proton <- function(x) {
  
  HX.SOBP.dose <- function(x) {
    
    dose.fun    <- function(idx, z) approx(DDD.list.sub[[idx]]$z.cm, DDD.list.sub[[idx]]$D.Gy, xout=z)$y
    
    for(n in 1:length(jj)) {
      if(n==1) { sum <- x[n] * dose.fun(idx=n, z=depth.seq) }
      else     { sum <- sum + x[n] * dose.fun(idx=n, z=depth.seq) }
    }
    
    return(sum)
    
  }
  
  dev            <- sum( abs(HX.SOBP.dose(abs(x))[min.depth.step:max.depth.step] - mean(HX.SOBP.dose(abs(x))[min.depth.step:max.depth.step])) ) * 1E7
  
  return(dev)
  
}