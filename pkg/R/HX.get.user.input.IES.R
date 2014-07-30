HX.get.user.input.IES <- function( n.IES,
                                   user.input,
                                   default.input,
                                   df.particles,
                                   df.fields,
                                   df.libc.HIT,
                                   input.con){

  # choose energy
  user.input$IES[[n.IES]]$energy.value.MeV.u <- HX.prompt.user( prompt        = "Select an (approximate) energy value in MeV/u",
                                                                default       = default.input$IES[[1]]$energy.value.MeV.u,
                                                                input.con     = input.con)
  
  user.input$IES[[n.IES]]$energy.value.MeV.u <- as.numeric(user.input$IES[[n.IES]]$energy.value.MeV.u)
  
  # Check for available energies in the libC
  if(!user.input$IES[[n.IES]]$energy.value.MeV.u%in%df.libc.HIT$energy.MeV.u){ 
    # if user value is not in libC suggest closest values
    user.input$IES[[n.IES]]$chosen.idx <- as.numeric( HX.prompt.user( prompt        = "Choose an energy value in MeV/u",
                                                                      choices       = HX.get.available.energies(user.input$IES[[n.IES]]$energy.value.MeV.u, df.libc.HIT),
                                                                      default       = default.input$IES[[1]]$chosen.idx,
                                                                      input.con     = input.con))
    
    user.input$IES[[n.IES]]$energy.value.MeV.u       <- HX.get.available.energies(user.input$IES[[n.IES]]$energy.value.MeV.u, df.libc.HIT)[user.input$IES[[n.IES]]$chosen.idx]
  }
  
  ### TODO ###
  ### MAKE FOLLOWING AVAILABLE TO NON-INTERACTIVE USE ALSO ###
  
  # Choose corresponding focus value
  energy.idx         <- which(df.libc.HIT$energy.MeV.u == user.input$IES[[n.IES]]$energy.value.MeV.u)
  foci.mm            <- unlist(df.libc.HIT[energy.idx,3:8])			# focus values from LIBC HIT table
  user.input$IES[[n.IES]]$chosen.foc.idx  <- as.numeric( HX.prompt.user( prompt        = "Choose focus (FWHM, mm) - PLEASE NOTE: F5 & F6 may be possible but are discouraged to use!",
                                                                         choices       = c(foci.mm, "auto"),
                                                                         default       = default.input$IES[[1]]$chosen.foc.idx,
                                                                         input.con     = input.con))
  if(user.input$IES[[n.IES]]$chosen.foc.idx == length(foci.mm)+1){
    user.input$IES[[n.IES]]$focus.FWHM.mm      <- 0.0
  }else{
    user.input$IES[[n.IES]]$focus.FWHM.mm      <- foci.mm[user.input$IES[[n.IES]]$chosen.foc.idx]
  }
  
  # Get intensity level
  user.input$IES[[n.IES]]$spill.intensity <- 0
  print(user.input$basic$intensity)
  if(user.input$basic$intensity > 1){ # 0 - not applicable, 1 - let system choose, i.e. omit node in plan
    user.input$IES[[n.IES]]$spill.intensity <- unlist(df.libc.HIT[energy.idx,9:23])[user.input$basic$intensity+1]			# intensity values from LIBC HIT table
  }
  
   # Choose field shape
   user.input$IES[[n.IES]]$field.shape.idx <- as.numeric( HX.prompt.user( prompt        = "What is the desired field shape?",
                                                                         choices       = as.character(df.fields$shape),
                                                                         default       = default.input$IES[[1]]$field.shape.idx,
                                                                         input.con     = input.con))

   
   # Choose size of homogenous field
   user.input$IES[[n.IES]]$field.size.mm <- as.numeric(HX.prompt.user( prompt        = "Please give side length of square homogenous area in irradiation field in mm",
                                                                       default       = default.input$IES[[1]]$field.size.mm,
                                                                       input.con     = input.con))

   #TODO: branch for circle (r.min oder so)
   if(df.fields$shape[user.input$IES[[n.IES]]$field.shape.idx] == "circle"){
     user.input$IES[[n.IES]]$r.min.mm <- as.numeric(HX.prompt.user( prompt        = "Please give inner radius of irradiation field in mm",
                                                                    default       = default.input$IES[[1]]$r.min.mm,
                                                                    input.con     = input.con))
   }
   
   # Choose fluence or dose
   user.input$IES[[n.IES]]$fluence.cm2.or.dose.Gy <- as.numeric(HX.prompt.user( prompt        = "Fluence (if positive number, in 1/cm2) or dose (negative, in Gy) for this field",
                                                                                default       = default.input$IES[[1]]$fluence.cm2.or.dose.Gy,
                                                                                input.con     = input.con))
  
  
  return( user.input )
} 
