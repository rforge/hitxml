rm(list =ls())
library(lattice)

files  <- list.files(pattern = "*_dEdx", path = ".", full.names = TRUE)

ll     <- lapply(1:length(files),
                 function(i,f){ tmp <- read.table(f[i],
                                                  header = TRUE,
                                                  col.names = c("name", 
                                                                "Z", 
                                                                "E.MeV.u", 
                                                                "dE.dx.keV.um"))[,2:4]
                                # keep proton and helium data only once
                                if(i != 1){
                                  ii  <- tmp$Z%in%c(1,2)
                                  tmp <- tmp[!ii,]
                                }
                                tmp},
                 files)

df      <- do.call("rbind", ll)
save(df, file = "dEdX_FLUKA_WATER.Rda")

xyplot( dE.dx.keV.um ~ E.MeV.u,
        df,
        type   = "l",
        grid   = TRUE,
        groups = Z,
        auto.key = list(title = "Z", 
                        space = "right", 
                        lines = TRUE, 
                        points = FALSE),
        scales = list(log = 10),
        ylab   = "stopping power / (keV/um)",
        xlab   = "kinetic energy / (MeV/u)",
        main   = "sg74001 dEdx customized water (FLUKA)")