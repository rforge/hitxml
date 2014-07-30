HX.save.defaults <- function(input){
  
  # A bit cumbersome, to keep backward compatibility
  df.defaults <-  data.frame( variable = c( "date",
                                            "name.exp.series",
                                            "name.exp.run",
                                            "chosen.particle",
                                            "chosen.rifi.idx",
                                            "intensity",
                                            "resolution.mm",
                                            "energy.value.MeV.u",
                                            "chosen.idx",
                                            "chosen.foc.idx",
                                            "field.shape.idx",
                                            "fluence.cm2.or.dose.Gy",
                                            "field.size.mm",
                                            "r.min.mm"),
                              value    = c( input$basic$date,
                                            input$basic$name.exp.series,
                                            input$basic$name.exp.run,
                                            input$basic$chosen.particle,
                                            input$basic$chosen.rifi.idx,
                                            input$basic$intensity,
                                            input$basic$resolution.mm,
                                            input$IES[[1]]$energy.value.MeV.u,
                                            input$IES[[1]]$chosen.idx,
                                            input$IES[[1]]$chosen.foc.idx,
                                            input$IES[[1]]$field.shape.idx,
                                            input$IES[[1]]$fluence.cm2.or.dose.Gy,
                                            input$IES[[1]]$field.size.mm,
                                            input$IES[[1]]$r.min.mm),
                              stringsAsFactors = FALSE)
        
    write.table( df.defaults,
                 system.file( "exec",
                              "HIT_XML.defaults",
                              package = "HITXML"),
                 quote            = FALSE,
                 row.names        = FALSE,
                 sep              = ";")

}