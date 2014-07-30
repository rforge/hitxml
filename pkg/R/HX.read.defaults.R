HX.read.defaults <- function(){
  
  df.defaults <- read.table( system.file( "exec",
                                         "HIT_XML.defaults",
                                         package = "HITXML"),
                            header  = TRUE,
                            sep              = ";",
                            stringsAsFactors = FALSE)

  var.value <- function(df, var.name){
    df$value[df$variable == var.name]
  }
  
  return( list( basic = list( date             = var.value(df.defaults, "date"),
                              name.exp.series  = var.value(df.defaults, "name.exp.series"),
                              name.exp.run     = var.value(df.defaults, "name.exp.run"),
                              chosen.particle  = var.value(df.defaults, "chosen.particle"),        
                              chosen.rifi.idx  = var.value(df.defaults, "chosen.rifi.idx"),
                              intensity        = var.value(df.defaults, "intensity"),
                              resolution.mm    = var.value(df.defaults, "resolution.mm")),
                                                          
                IES   = list( list( energy.value.MeV.u     = var.value(df.defaults, "energy.value.MeV.u"),
                                    chosen.idx             = var.value(df.defaults, "chosen.idx"),
                                    chosen.foc.idx         = var.value(df.defaults, "chosen.foc.idx"),
                                    field.shape.idx        = var.value(df.defaults, "field.shape.idx"),      
                                    fluence.cm2.or.dose.Gy = var.value(df.defaults, "fluence.cm2.or.dose.Gy"),
                                    field.size.mm          = var.value(df.defaults, "field.size.mm"),
                                    r.min.mm               = var.value(df.defaults, "r.min.mm")))))     
}