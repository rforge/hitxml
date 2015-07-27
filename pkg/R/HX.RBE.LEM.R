# Computes the relative biological effectiveness (RBE) according to
# Krämer et al. (2000) using statistical sampling
#
# Arguments
#  RBE.data       dataRBE object
#  N.event        number of samples for averaging
#  D.Gy           dose chosen for the target volume
#  spc            spc-formatted data file with energy-fluence in n depth steps
#                 (AT.SPC.read(file.name, endian, flavour))
#  df.spect       spc data of one depth step
#
# Return value
#  A matrix with RBE,
#                biological dose 
#                and absorbed dose

HX.RBE.LEM <- function (RBE.data, Spectrum.data, dose.Gy, N.events = 100, N.runs = 10, write.output = FALSE){
  
  ##################################################
  # calculate the number of hits on a cell nucleus #
  ##################################################
  A.nucl.um2      <- pi * (r.nucleus.um(RBE.data))^2 
  A.nucl.cm2      <- A.nucl.um2 / (10000^2)
  
  # relative fluence (particle per primary particle)
  fluence.spectrum <- particles.per.primary(Spectrum.data) 
  dose.spectrum.Gy <- dose.per.primary(Spectrum.data, "PSTAR")
  # fluence factor to get dose set
  fluence.factor   <- dose.Gy / dose.spectrum.Gy
  
  
  # mean number of hits on a cell nucleus (N.Hit.average)
  N.hit.avg <- A.nucl.cm2 * fluence.spectrum * fluence.factor
  
  # N.Hit has to be sampled at random from the Poisson distribution
  N.hit <- rpois(N.events, N.hit.avg)
   
  # Pre-compute stopping powers
  S.MeV.cm2.g  <- Mass.Stopping.Power.MeV.cm2.g(Spectrum.data, "PSTAR")
  ###############################################################
  # create a set of particles (particle type T(k), energy E(k)) #
  ###############################################################
  
  D.abs.Nhit.Gy <- numeric(N.events) # initialization (absorbed dose)
  N.lethal.Nhit <- numeric(N.events) # initialization (effective damage)
  
  # index
  idx          <- 1:nrow(Spectrum.data@spectrum)
  
  results      <- NULL
  # Loop over all runs
  for(j in 1:N.runs){
    # Loop over all events
    # j <- 1
    for (i in 1:N.events){
      # i <- 1
      
      # in case no particle hits a cell nucleus, everything is 0
      if (N.hit[i] == 0) {
        D.abs.Nhit.Gy[i] <- 0
        N.lethal.Nhit[i] <- 0
      }else{
        
        # sampled n = N.Hit times to obtain a set (T(k), E(k))
        particles.idx <- sample(x       = idx,
                                size    = N.hit[i],
                                replace = TRUE,
                                prob    = Spectrum.data@spectrum$N.per.primary)
        
        ##############################
        # new table for sampled data #
        ##############################
        
        particle.no  <- Spectrum.data@spectrum$particle.no[particles.idx]
        E.MeV.u      <- Spectrum.data@spectrum$E.mid.MeV.u[particles.idx]
        LET          <- S.MeV.cm2.g[particles.idx]
        
        
        #################################################################
        # compute absorbed dose (D_abs) and effective damage (N_lethal) #
        #################################################################
        
        # (1) absorbed dose 
        
        c.1.cm2            <- 1.60217657e-10 / A.nucl.cm2
        d.Gy               <- LET * c.1.cm2
        D.abs.Gy           <- cumsum(d.Gy) # returns a vector whose elements are the cumulative sums
        
        
        # (2) effective damage
        
        # D.abs(k-1)
        D.k1               <- D.abs.Gy - d.Gy
        
        jj                 <- D.k1 < D.cut.Gy(RBE.data)
        
        # alpha.ion for all particles.idx
        cur.alpha.ion      <- alpha.ion(RBE.data, 
                                        AT.particle.name.from.particle.no(particle.no), 
                                        E.MeV.u)
        
        # dose dependent slopes of the heavy ion effect curve
        s.1.Gy      <- cur.alpha.ion + (s.max(RBE.data) - cur.alpha.ion) * D.k1 / D.cut.Gy(RBE.data)
        s.1.Gy[!jj] <- s.max(RBE.data)
        
        # effective damage
        n        <- s.1.Gy * LET * c.1.cm2
        N.lethal <- cumsum(n)
        
        
        # D.abs and N.lethal for N.hit
        D.abs.Nhit.Gy[i] <- D.abs.Gy[N.hit[i]]
        N.lethal.Nhit[i] <- N.lethal[N.hit[i]]
        
        # Output progress
        if(write.output){
          cat("Done event no. ", i, " in run ", j, 
              "- hits: ", N.hit[i], 
              " | D.abs.Gy: ", D.abs.Gy[N.hit[i]],
              " | lethal hits: ", N.lethal[N.hit[i]],
              "\n")
        }
        
      }
      
    } # N.event
    
    
    #average absorbed dose
    D.abs.average.Gy <- sum (D.abs.Nhit.Gy) / N.events
    
    #average survival
    s.average.1.Gy <- sum( exp( -1.0*N.lethal.Nhit)) / N.events
    
    #average damage
    N.lethal.average <- -log(s.average.1.Gy)
    
    
    ###################
    # biological dose #
    ###################
    
    # compute D.biol by inversion of the linear-quadratic equation for the cellular survival 
    
    u <- alpha.X(RBE.data) / (2 * beta.X(RBE.data))
    
    # case differentiation: D < D.cut & D >= D.cut
    D1.Gy <- - u + sqrt(u^2 + N.lethal.average / beta.X(RBE.data))
    D2.Gy <- (N.lethal.average - alpha.X(RBE.data) * D.cut.Gy(RBE.data) - beta.X(RBE.data) * D.cut.Gy(RBE.data)^2) / s.max(RBE.data) + D.cut.Gy(RBE.data)
    
    if (D1.Gy < D.cut.Gy(RBE.data)){
      D.biol.Gy <- D1.Gy
    }else{
      D.biol.Gy <- D2.Gy}
    
    
    ###############
    # compute RBE #
    ###############
    
    LEM.RBE <- D.biol.Gy / D.abs.average.Gy
    
    run.results <- matrix(ncol=3, nrow=1, data = c(LEM.RBE, D.biol.Gy, D.abs.average.Gy))
    results     <- rbind(results, run.results)
  }# N.runs
  
  
  return(c( RBE             = mean(results[,1]),
            u.rel.RBE       = sd(results[,1]) / sqrt(N.runs) / mean(results[,1]),
            D.biol.Gy       = mean(results[,2]),
            u.rel.D.biol.Gy = sd(results[,2]) / sqrt(N.runs) / mean(results[,2]),
            D.Gy            = mean(results[,3]),
            u.rel.D.Gy      = sd(results[,3]) / sqrt(N.runs) / mean(results[,3])))
}