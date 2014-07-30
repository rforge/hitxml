HX.correct.user.input <- function(user.input, n.IES, df.libc.HIT){
  # Energy
  if(!user.input$IES[[n.IES]]$energy.value.MeV.u%in%df.libc.HIT$energy.MeV.u){ 
    # if user value is not in libC suggest closest values
    user.input$IES[[n.IES]]$chosen.idx           <- 1
    user.input$IES[[n.IES]]$energy.value.MeV.u   <- HX.get.available.energies(user.input$IES[[n.IES]]$energy.value.MeV.u, df.libc.HIT)[1]
  }
  
  energy.idx         <- which(df.libc.HIT$energy.MeV.u == user.input$IES[[n.IES]]$energy.value.MeV.u)
  foci.mm            <- unlist(df.libc.HIT[energy.idx,3:8])  		# focus values from LIBC HIT table
  user.input$IES[[n.IES]]$chosen.foc.idx  <- as.numeric(default.input$IES[[1]]$chosen.foc.idx)
  
  if(user.input$IES[[n.IES]]$chosen.foc.idx == length(foci.mm)+1){
    user.input$IES[[n.IES]]$focus.FWHM.mm      <- 0.0
  }else{
    user.input$IES[[n.IES]]$focus.FWHM.mm      <- foci.mm[user.input$IES[[n.IES]]$chosen.foc.idx]
  }
  
  # Get intensity level
  user.input$IES[[n.IES]]$spill.intensity <- 0
  if(user.input$basic$intensity > 1){ # 0 - not applicable, 1 - let system choose, i.e. omit node in plan
    user.input$IES[[n.IES]]$spill.intensity <- unlist(df.libc.HIT[energy.idx,9:23])[user.input$basic$intensity+1]			# intensity values from LIBC HIT table
  }
  
  # Choose field shape
  user.input$IES[[n.IES]]$field.shape.idx <- as.numeric(default.input$IES[[1]]$field.shape.idx)
  
  # Choose size of homogenous field
  user.input$IES[[n.IES]]$field.size.mm   <- as.numeric(default.input$IES[[1]]$field.size.mm)
  
  #TODO: branch for circle (r.min oder so)
#  if(field.shape[field.shape.idx] == "circle"){
#    user.input$IES[[n.IES]]$r.min.mm <- as.numeric(default.input$IES[[1]]$r.min.mm)
#  }
  
  # Choose fluence or dose
  user.input$IES[[n.IES]]$fluence.cm2.or.dose.Gy <- as.numeric(default.input$IES[[1]]$fluence.cm2.or.dose.Gy)
  
  return(user.input)  
}