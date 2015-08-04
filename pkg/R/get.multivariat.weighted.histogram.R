#' @title Multivariate histogram from weighted variable
#' 
#' @description Creates histogram for a continous, grouped variable that can be in 
#' addition be weighted.
#' 
#' @param x continous variable (vector)
#' @param i variable to group i (vector of same length)
#' @param w optional weighting variable
#' @param x.min minimal value (left side for lowest bin)
#' @param x.max maximal value (right side for highest bin)
#' @param n.x.bin number of bins
#' @param log if true, binning will be done on a logarithmic scale
#' 
#' @return data frame with bin mids, bin widths and density
#' 
#' @details Bins are left sided. Values outside [x.min, x.max) are disregarded with a warning
#' 
#' @author Steffen Greilich
get.multivariate.weighted.histogram <- function(x, i, w, x.min, x.max, n.x.bins, log = FALSE){

  i.set    <- sort(unique(i))
  i.max    <- max(i.set)
  n.i      <- length(i.set)
  
  if(missing(w)){
    w <- rep(1.0, length(x))
  }
  
  # Set-up histogram
  if(log){
    breaks   <- 10^seq( from       = log10(x.min), 
                        to         = log10(x.max), 
                        length.out = n.x.bins+1)
  }else{
    breaks   <- seq( from       = x.min, 
                     to         = x.max, 
                     length.out = n.x.bins+1)
  }
  
  df.hist  <- data.frame(x.low     = rep(breaks[-length(breaks)], n.i),
                         x.high    = rep(breaks[-1],              n.i),
                         i         = sort(rep(i.set, n.x.bins)),
                         frequency = numeric((length(breaks)-1) * n.i))
  
  for (cur.i in i.set){
    # cur.i <- 1
    ii  <- i == cur.i
    jj  <- df.hist$i == cur.i
    for(j in 1:sum(jj)){
      # k <- 1
      kk <- x[ii] >= df.hist$x.low[jj][j] & x[ii] < df.hist$x.high[jj][j]
      df.hist$frequency[jj][j] <- sum(w[ii][kk])
    }
    cat("Done i =", cur.i, "\n")
  }
  
  if(log){
    df.hist$x             <- sqrt(df.hist$x.high * df.hist$x.low)
  }else{
    df.hist$x             <- 0.5*(df.hist$x.high + df.hist$x.low)
  }
  df.hist$dx            <- df.hist$x.high - df.hist$x.low
  df.hist               <- df.hist[,!names(df.hist)%in%c("x.low", "x.high")]
  
  df.hist$density       <- df.hist$frequency / df.hist$dx
  
  return(df.hist)
}