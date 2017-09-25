plot.SOBP <- function(plot.ddds, plot.depths.g.cm2, plot.weights, add.comment, start.depth.cm = 0, end.depth.cm = 0){
  plot.no.IES             <- length(plot.ddds@beam.energies.MeV.u)
  
  D.Gy                    <- get.dose.Gy.from.set( DDD.set      = plot.ddds, 
                                                   depths.g.cm2 = plot.depths.g.cm2, 
                                                   weights      = plot.weights)

  dd                      <- unlist( lapply( 1:plot.no.IES,
                                             function(i){
                                                get.dose.Gy(get.ddd(plot.ddds, 
                                                                    ddds.sub@beam.energies.MeV.u[i]), 
                                                            plot.depths.g.cm2) * plot.weights[i]}))
  
  df <- data.frame( depths.g.cm2 = rep(plot.depths.g.cm2, plot.no.IES+1),
                    D.Gy         = c(D.Gy, dd),
                    which        = c(rep("SOBP", length(plot.depths.g.cm2)),
                                     sort(rep(1:plot.no.IES, length(plot.depths.g.cm2)))))
  xyplot(D.Gy ~ depths.g.cm2,
         df,
         grid = TRUE,
         type = "l",
         col  = c(rep("grey", plot.no.IES), "blue"),
         alpha  = c(rep(0.8, plot.no.IES), 1),
         groups = which,
         xlab = list("distal depth /cm", cex=1.5),
         ylab = list("total dose / Gy", cex=1.5),
         scale = list(cex = 1.25),
         panel = function(...){
                  panel.xyplot(...)
                  panel.abline(v = start.depth.cm, lty = 2)
                  panel.abline(v = end.depth.cm, lty = 2)
         },
         main = list(paste0("SOBP (single field, ",
                            unique(ddds.sub@projectiles),
                            ") consisting of ", 
                            plot.no.IES, 
                            " IESs ",
                            add.comment), 
                     cex=1.5),
         sub  = "NB: HIT isocenter is at 0.289 g/cm2 depth!")
}