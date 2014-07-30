

# Check integration
    if(F){
        PBR                <- data.frame( x.mm       = c(0,15),
                                           y.mm       = c(0,15),
                                           focusX.mm  = c(10,5),
                                           focusY.mm  = c(10,5),
                                           particles  = c(0.67, 0.33),
                                           expid      = 'optim')
    m <- PBR.compute.field(PBR           = PBR, 
                           expid         = 'optim', 
                           resolution.mm = 0.1)
    levelplot(m[,3] ~ m[,1]*m[,2])
    sum(m[,3]*0.1*0.1)
    }
