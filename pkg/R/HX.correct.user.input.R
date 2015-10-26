HX.correct.user.input <- function(user.input, n.IES, df.libc.HIT){
  # Energy
  if(is.null(user.input$IES[[n.IES]]$chosen.idx)){
    user.input$IES[[n.IES]]$energy.value.MeV.u   <- HX.get.available.energies(user.input$IES[[n.IES]]$energy.value.MeV.u, df.libc.HIT)[1]
    user.input$IES[[n.IES]]$chosen.idx           <- df.libc.HIT$energy.no[match(user.input$IES[[n.IES]]$energy.value.MeV.u, df.libc.HIT$energy.MeV.u)]
  }else{
    user.input$IES[[n.IES]]$energy.value.MeV.u   <- df.libc.HIT$energy.value.MeV.u[match(user.input$IES[[n.IES]]$chosen.idx, df.libc.HIT$energy.no)]    
  }
    

  # Foci
  idx                <- match(user.input$IES[[n.IES]]$chosen.idx, df.libc.HIT$energy.no)
  foci.mm            <- unlist(df.libc.HIT[idx,3:8])  		# focus values from LIBC HIT table
  if(is.null(user.input$IES[[n.IES]]$chosen.foc.idx)){
    dd                                           <- (foci.mm - user.input$IES[[n.IES]]$focus.FWHM.mm)^2
    user.input$IES[[n.IES]]$chosen.foc.idx       <- which(min(dd) == dd)
    user.input$IES[[n.IES]]$focus.FWHM.mm        <- foci.mm[user.input$IES[[n.IES]]$chosen.foc.idx]
  }else{
    user.input$IES[[n.IES]]$focus.FWHM.mm        <- foci.mm[user.input$IES[[n.IES]]$chosen.foc.idx]    
  }
  
  # Get intensity level
  user.input$IES[[n.IES]]$spill.intensity <- 0
  if(user.input$basic$intensity > 1){ # 0 - not applicable, 1 - let system choose, i.e. omit node in plan
    user.input$IES[[n.IES]]$spill.intensity <- unlist(df.libc.HIT[energy.idx,9:23])[user.input$basic$intensity+1]			# intensity values from LIBC HIT table
  }
  
  # Choose field shape
#  user.input$IES[[n.IES]]$field.shape.idx <- as.numeric(default.input$IES[[1]]$field.shape.idx)
  
  # Choose size of homogenous field
#  user.input$IES[[n.IES]]$field.size.mm   <- as.numeric(default.input$IES[[1]]$field.size.mm)
  
  #TODO: branch for circle (r.min oder so)
#  if(field.shape[field.shape.idx] == "circle"){
#    user.input$IES[[n.IES]]$r.min.mm <- as.numeric(default.input$IES[[1]]$r.min.mm)
#  }
  
  # Choose fluence or dose
#  user.input$IES[[n.IES]]$fluence.cm2.or.dose.Gy <- as.numeric(default.input$IES[[1]]$fluence.cm2.or.dose.Gy)
  
  return(user.input)  
}