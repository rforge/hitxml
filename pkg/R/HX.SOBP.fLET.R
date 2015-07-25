HX.SOBP.fLET <- function(x) {
  
  fLET.fun    <- function(idx, z) approx(df.total[[idx]]$depth.g.cm2, df.total[[idx]]$fLET.keV.um, xout=z)$y
  
  for(n in 1:length(jj)) {
    if(n==1) { sum <- x[n] * fLET.fun(idx=n, z=depth.seq) }
    else     { sum <- sum + x[n] * fLET.fun(idx=n, z=depth.seq) }
  }
  
  return(sum)
  
}
