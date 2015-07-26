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

HX.RBE.LEM <- function (RBE.data, N.event, N.runs, D.Gy, spc, df.spect, write.output = FALSE){
  
  ##################################################
  # calculate the number of hits on a cell nucleus #
  ##################################################
  A.nucl.um2      <- pi * (r.nucleus.um(RBE.data))^2 
  A.nucl.cm2      <- A.nucl.um2 / (10000^2)
  
  # relative fluence (particle per primary particle)
  rel.fluence.sum <- sum(df.spect$N.per.primary) 
  
  # calculation of the phys. dose at the position of the Bragg peak
  df.spect.1  <- spc$spc[spc$spc$depth.step == 1,]
  #df.spect.45 <- spc$spc[spc$spc$depth.step == 45,] # depthstep 45: position of the Bragg peak
  # [can be checked with:
  #  which(Results$D.abs.Gy == max(Results$D.abs.Gy))]
  #df.spect.45$N.per.primary <- df.spect.45$N.per.primary / sum(df.spect.1$N.per.primary) # normalization
  
  
  #dose.depth45.Gy <- sum(AT.dose.Gy.from.fluence.cm2(E.MeV.u = df.spect.45$E.mid.MeV.u,       
  #particle.no = df.spect.45$particle.no,   
  #fluence.cm2 = df.spect.45$N.per.primary,
  #material.no = 1,
  #stopping.power.source.no = 0)$dose.Gy)
  # AT.dose.Gy.from.fluence.cm2: returns dose in Gy for each given particle
  
  # primary fluence for the given dose
  #prim.fluence <- D.Gy / (dose.depth45.Gy)
  
  #effective fluence
  eff.fluence.1.cm2 <- rel.fluence.sum  * 217647.3  # 4Gy: F= 43529457, 0.02Gy: F = 217647.3
  
  
  # mean number of hits on a cell nucleus (N.Hit.average)
  N.hit.avg <- A.nucl.cm2 * eff.fluence.1.cm2
  
  # N.Hit has to be sampled at random from the Poisson distribution
  N.hit <- rpois(N.event, N.hit.avg)
  
  
  #######################################################
  # convert the fluence into a probability distribution #
  #######################################################
  
  df.spect$norm.fluence <- df.spect$N.per.primary / sum(df.spect$N.per.primary)
  
  
  #################################
  # compute (mass) stopping power #
  #################################
  
  #df.spect$S.MeV.cm2.g <- AT.Stopping.Power.MeV.cm2.g(stopping.power.source.no = 0,
  #                                                    E.MeV.u                  = df.spect$E.mid.MeV.u, 
  #                                                    particle.no              = df.spect$particle.no, 
  #                                                    material.no              = 1)$Stopping.Power.MeV.cm2.g
  # AT.Stopping.Power.MeV.cm2.g: main access method to stopping power data
  
  ###############################################################
  # create a set of particles (particle type T(k), energy E(k)) #
  ###############################################################
  
  D.abs.Nhit.Gy <- numeric(N.event) # initialization (absorbed dose)
  N.lethal.Nhit <- numeric(N.event) # initialization (effective damage)
  
  # index
  df.spect$idx <- 1:nrow(df.spect)
  
  results      <- NULL
  # Loop over all runs
  for(j in 1:N.runs){
    # Loop over all events
    for (i in 1:N.event){
      
      # first case: no particle hits a cell nucleus
      if (N.hit[i] == 0) {
        D.abs.Nhit.Gy[i] <- 0
        N.lethal.Nhit[i] <- 0
      }
      
      # second case: at least one particle hits a cell nucleus
      else{
        
        # sampled n = N.Hit times to obtain a set (T(k), E(k))
        particles.idx <- sample(x = df.spect$idx,
                                size = N.hit[i],
                                replace = TRUE,
                                prob = df.spect$norm.fluence)
        
        ##############################
        # new table for sampled data #
        ##############################
        
        df.sample <- data.frame( #E.MeV.u      = df.spect$E.mid.MeV.u[particles.idx],
          particle.no  = df.spect$particle.no[particles.idx],
          S.MeV.cm2.g  = df.spect$S.MeV.cm2.g[particles.idx])
        
        
        #################################################################
        # compute absorbed dose (D_abs) and effective damage (N_lethal) #
        #################################################################
        
        # (1) absorbed dose 
        
        c.1.cm2  <- 1.6 * 10^(-10) * (A.nucl.cm2)^(-1)
        df.sample$d.Gy     <- df.sample$S.MeV.cm2.g * c.1.cm2
        df.sample$D.abs.Gy <- cumsum(df.sample$d.Gy) # returns a vector whose elements are the cumulative sums
        
        
        # (2) effective damage
        
        # D.abs(k-1)
        df.sample$D.k1 <- df.sample$D.abs.Gy - df.sample$d.Gy
        
        # boolean vector
        df.sample$jj <- df.sample$D.k1 < get.D.cut.Gy(RBE.data, projectile = spc$projectile)
        
        # alpha.ion for all particles.idx
        alpha.ion.1.Gy <- numeric(length(particles.idx))
        alpha.fun <- function(x) {
          get.alpha.ion.1.Gy(RBE.data = RBE.data,
                             projectile = AT.particle.name.from.particle.no(df.spect$particle.no[particles.idx[x]]),
                             S.MeV.cm2.g = df.spect$S.MeV.cm2.g[particles.idx[x]])
        }
        alpha.ion.1.Gy <- as.numeric(lapply(c(1:length(particles.idx)), alpha.fun))
        
        # dose dependent slopes of the heavy ion effect curve
        df.sample$s.1.Gy = alpha.ion.1.Gy + (get.s.max.1.Gy(RBE.data, projectile= spc$projectile) - alpha.ion.1.Gy) *
          df.sample$D.k1 / get.D.cut.Gy(RBE.data, projectile = spc$projectile)
        df.sample$s.1.Gy[!df.sample$jj] <- get.s.max.1.Gy(RBE.data, projectile = spc$projectile)
        
        # effective damage
        df.sample$n        <- df.sample$s.1.Gy * df.sample$S.MeV.cm2.g * c.1.cm2
        df.sample$N.lethal <- cumsum(df.sample$n)
        
        
        # D.abs and N.lethal for N.hit
        D.abs.Nhit.Gy[i] <- df.sample$D.abs.Gy[N.hit[i]]
        N.lethal.Nhit[i] <- df.sample$N.lethal[N.hit[i]]
        
        # Output progress
        if(write.output){
          cat("Done event no. ", i, " in run ", j, 
              "- hits: ", N.hit[i], 
              " | D.abs.Gy: ", df.sample$D.abs.Gy[N.hit[i]],
              " | lethal hits: ", df.sample$N.lethal[N.hit[i]],
              "\n")
        }
        
      }
      
    } # N.event
    
    
    ######################################
    # new table for the different events #
    ######################################
    
    df.events <- data.frame( D.abs.Nhit.Gy = D.abs.Nhit.Gy,
                             N.lethal.Nhit = N.lethal.Nhit)
    
    
    ############
    # averages #
    ############
    
    #average absorbed dose
    D.abs.average.Gy <- sum (df.events$D.abs.Nhit.Gy) / N.event
    
    #average survival
    s.average.1.Gy <- sum( exp( - df.events$N.lethal.Nhit)) / N.event
    
    #average damage
    N.lethal.average <- -log(s.average.1.Gy)
    
    
    ###################
    # biological dose #
    ###################
    
    # compute D.biol by inversion of the linear-quadratic equation for the cellular survival 
    
    u <- get.alpha.1.Gy(RBE.data, projectile= spc$projectile) /
      (2 * get.beta.1.Gy2(RBE.data, projectile= spc$projectile))
    
    # case differentiation: D < D.cut & D >= D.cut
    D1.Gy <- - u + sqrt(u^2 + N.lethal.average / get.beta.1.Gy2(RBE.data, projectile= spc$projectile))
    D2.Gy <- (N.lethal.average -
                get.alpha.1.Gy(RBE.data, projectile= spc$projectile) *
                get.D.cut.Gy(RBE.data, projectile= spc$projectile) -
                get.beta.1.Gy2(RBE.data, projectile= spc$projectile) *
                (get.D.cut.Gy(RBE.data, projectile= spc$projectile))^2) /
      get.s.max.1.Gy(RBE.data, projectile= spc$projectile) +
      get.D.cut.Gy(RBE.data, projectile= spc$projectile)
    
    # boolean vector
    ii <- D1.Gy < get.D.cut.Gy(RBE.data, projectile= spc$projectile)
    
    if (ii) {D.biol.Gy <- D1.Gy} else{D.biol.Gy <- D2.Gy}
    
    
    ###############
    # compute RBE #
    ###############
    
    LEM.RBE <- D.biol.Gy / D.abs.average.Gy
    
    run.results <- matrix(ncol=3, nrow=1, data = c(LEM.RBE, D.biol.Gy, D.abs.average.Gy))
    results <- rbind(results, run.results)
  }# N.runs
  
  df <- data.frame( RBE         = mean(results[,1]),
                    sd.RBE      = sd(results[,1]),
                    u.RBE       = sd(results[,1]) / sqrt(N.runs),
                    D.biol.Gy   = mean(results[,2]),
                    sd.D.biol.Gy = sd(results[,2]),
                    u.D.biol.Gy = sd(results[,2]) / sqrt(N.runs),
                    D.Gy        = mean(results[,3]),
                    sd.D.Gy     = sd(results[,3]),
                    u.D.Gy      = sd(results[,3]) / sqrt(N.runs))
  return(df)
}