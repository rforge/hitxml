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
  
  # Get parameters from RBE.data
  cur.alpha.X     <- alpha.X(RBE.data)
  cur.beta.X      <- beta.X(RBE.data)
  cur.D.cut.Gy    <- D.cut.Gy(RBE.data)
  cur.s.max       <- s.max(RBE.data)
  
  # calculate the number of hits on a cell nucleus #
  A.nucl.um2      <- pi * (r.nucleus.um(RBE.data))^2 
  A.nucl.cm2      <- A.nucl.um2 / (10000^2)
  c.1.cm2         <- 1.60217657e-10 / A.nucl.cm2
  
  # relative fluence (particle per primary particle)
  fluence.spectrum <- particles.per.primary(Spectrum.data) 
  dose.spectrum.Gy <- dose.per.primary(Spectrum.data, "PSTAR")

  # fluence factor to get dose set
  fluence.factor   <- dose.Gy / dose.spectrum.Gy
  
  # mean number of hits on a cell nucleus (N.Hit.average)
  N.hit.avg <- A.nucl.cm2 * fluence.spectrum * fluence.factor
  
  # Pre-compute stopping powers
  S.MeV.cm2.g  <- Mass.Stopping.Power.MeV.cm2.g(Spectrum.data, "PSTAR")

  # index for faster particle sampling
  idx          <- 1:nrow(Spectrum.data@spectrum)
  
  results      <- NULL
  
  # Sampled number of hits, get indices
  N.hit         <- rpois(N.events * N.runs, N.hit.avg)
    
  i <- unlist( mapply( function(x,n){ rep(x,n)},
                       x = seq_along(N.hit),
                       n = N.hit))

  j <- unlist( mapply( function(x,n){ if(n != 0){
                                        rep(x,n)
                                        }else{
                                          numeric(0)
                                        }},
                       x = sort(rep(1:N.runs, N.events)),
                       n = N.hit))

  # sampled N.hit sets of (T(k), E(k))
  particles.idx <- sample(x       = idx,
                          size    = length(i),
                          replace = TRUE,
                          prob    = Spectrum.data@spectrum$N.per.primary)

  particle.no   <- Spectrum.data@spectrum$particle.no[particles.idx]
  E.MeV.u       <- Spectrum.data@spectrum$E.mid.MeV.u[particles.idx]
  LET           <- S.MeV.cm2.g[particles.idx]


  # Compute absorbed dose 
  d.Gy          <- LET * c.1.cm2
  D.abs.Gy      <- unlist( tapply(d.Gy, i, cumsum))
    

  # Compute effective damage
  D.k1          <- D.abs.Gy - d.Gy
  jj            <- D.k1 < cur.D.cut.Gy
        
  # alpha.ion for all particles.idx
  cur.alpha.ion <- alpha.ion( RBE.data, 
                                AT.particle.name.from.particle.no(particle.no), 
                                E.MeV.u)
        
  # dose dependent slopes of the heavy ion effect curve
  s.1.Gy        <-  cur.alpha.ion + (cur.s.max - cur.alpha.ion) * D.k1 / cur.D.cut.Gy
  s.1.Gy[!jj]   <- cur.s.max

  n             <- s.1.Gy * LET * c.1.cm2
        
  # D.abs and N.lethal for N.hit
  D.abs.Nhit.Gy <- tapply(D.abs.Gy, i, max)
  N.lethal.Nhit <- tapply(n, i, sum)
  j             <- tapply(j, i, unique)
  
  # average absorbed dose (for N.runs chunks to assess variance)
  D.abs.average.Gy <- tapply(D.abs.Nhit.Gy, j, function(x){sum(x)/N.events})
    
  # average survival
  s.average.1.Gy   <- tapply(N.lethal.Nhit, j, function(x){sum( exp( -1.0*x))/N.events})
    
  #average damage
  N.lethal.average <- -log(s.average.1.Gy)
    
  # compute D.biol by inversion of the linear-quadratic equation for the cellular survival 
  u <- cur.alpha.X / (2 * cur.beta.X)
    
  # case differentiation: D < D.cut or D >= D.cut?
  D1.Gy <- - u + sqrt(u^2 + N.lethal.average / cur.beta.X)
  D2.Gy <- (N.lethal.average - cur.alpha.X * cur.D.cut.Gy - cur.beta.X * cur.D.cut.Gy^2) / cur.s.max + cur.D.cut.Gy
    
  D.biol.Gy <- D1.Gy
  ii    <- D1.Gy >= cur.D.cut.Gy
  D.biol.Gy[ii] <- D2.Gy[ii]
    
  LEM.RBE <- D.biol.Gy / D.abs.average.Gy

  return(c( RBE             = mean(LEM.RBE),
            u.rel.RBE       = sd(LEM.RBE) / sqrt(N.runs) / mean(LEM.RBE),
            D.biol.Gy       = mean(D.biol.Gy),
            u.rel.D.biol.Gy = sd(D.biol.Gy) / sqrt(N.runs) / mean(D.biol.Gy),
            D.Gy            = mean(D.abs.average.Gy),
            u.rel.D.Gy      = sd(D.abs.average.Gy) / sqrt(N.runs) / mean(D.abs.average.Gy)))
}