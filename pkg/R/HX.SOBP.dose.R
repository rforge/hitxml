HX.SOBP.dose <- function(x) {
  
  dose.fun    <- function(idx, z) approx(DDD.list.sub[[idx]]$z.cm, DDD.list.sub[[idx]]$D.Gy, xout=z)$y
      
  for(n in 1:length(jj)) {
    if(n==1) { sum <- x[n] * dose.fun(idx=n, z=depth.seq) }
    else     { sum <- sum + x[n] * dose.fun(idx=n, z=depth.seq) }
  }
  
  return(sum)
  
}
