rm(list = ls())

library(HITXML)
library(libamtrack)

####################
# USER INPUT
####################

# input of the energy, the projectile and the material
tissue           <- "chordom02"
E.MeV.u          <- 210            ########## UNITS???
projectile       <- "12C"


##############################################################
# MODIFIED READING FUNCTION TO WORK ON SINGLE FILE,
# REMOVED COMPUTATION OF ANALYTIC STOPPING POWER (UNNEC.)
##############################################################
RBE.data <- dataRBE(tissue)




##########################################
# Compute stopping power
# EVERYTHING IS WATER...
##########################################
S.MeV.cm    <- AT.Mass.Stopping.Power( stopping.power.source = "ICRU",
                                       E.MeV.u               = E.MeV.u,
                                       particle.no           = AT.particle.no.from.particle.name(projectile),
                                       material.no           = AT.material.no.from.material.name("Water, Liquid"))$stopping.power.MeV.cm2.g

# compute alpha_ion, beta_ion
alpha_ion  <- alpha.1.Gy(RBE.data) * RBE.alpha(RBE.data, projectile, E.MeV.u)

###################################################
###################################################
###################################################
                          
                          
s_max      <- alpha.1.Gy(RBE.data) + 2 * beta.1.Gy2(RBE.data) * D.cut.Gy(RBE.data)

beta_ion   <- (s_max - alpha_ion)/(2 * D.cut.Gy(RBE.data))

# compute  RBE
A_nucl = pi * (r.nucleus.um(RBE.data))^2
ce      = 1.6 * 10^(-2) * (A_nucl)^(-1)

n = 1000    # number of samples (N_event)

N_lethal    <- c(0:n)
D_abs       <- c(0:n)
s           <- c(1:n)


for (k in 2:(n+1)){
  D_abs[k]    = D_abs[k-1] + S.MeV.cm *ce
}

for (k in 1:n){
  if (D_abs[k] < D.cut.Gy(RBE.data)){
    s[k] = alpha_ion + (s_max - alpha_ion) * D_abs[k] /D.cut.Gy(RBE.data)
    }
  else {s[k] = s_max}
}

for (k in 2:(n+1)){
  N_lethal[k] = N_lethal[k-1] + s[k-1] * S.MeV.cm *ce
}

 #average absorbed dose
D_abs_average = 0

for (i in 1:n){
  D_abs_average <- D_abs_average + D_abs[i+1]
}

D_abs_average <- D_abs_average /n


#average survival
s_average = 0

for (i in 1:n){
  s_average <- s_average + exp(-N_lethal[i+1])
}

s_average <- s_average /n


#average damage
N_lethal_average <- -log(s_average)


#compute D_biol

u = alpha.1.Gy(RBE.data) / (2 * beta.1.Gy2(RBE.data))

D1 = - u + sqrt(u^2 + N_lethal_average/beta.1.Gy2(RBE.data))
D2 = (N_lethal_average - alpha.1.Gy(RBE.data) * D.cut.Gy(RBE.data) - beta.1.Gy2(RBE.data) * (D.cut.Gy(RBE.data)^2)) / s_max + D.cut.Gy(RBE.data)

if (D1 < D.cut.Gy(RBE.data) & D2 < D.cut.Gy(RBE.data)) {D_biol = D1} else{D_biol = D2}

#compute RBE

RBE = D_biol/ D_abs_average
