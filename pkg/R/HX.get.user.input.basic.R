HX.get.user.input.basic <- function( user.input,
                                     default.input,
                                     df.particles, 
                                     df.rippleFilter,
                                     input.con){
  
  # Get experiment date
  if(default.input$basic$date == -1){
    default.input$basic$date      <- format(Sys.Date(), "%Y%m%d")
  }
  user.input$basic$date        <- HX.prompt.user( prompt    = "Please give date of experiment , YYYYMMDD",
                                                  default   = default.input$basic$date,
                                                  input.con = input.con)  
  
  # Get experiment series name
  user.input$basic$name.exp.series          <- HX.prompt.user( prompt    = "Please give experiment series name",
                                                               default   = default.input$basic$name.exp.series,
                                                               input.con = input.con)

  # Get experiment run name
  user.input$basic$name.exp.run          <- HX.prompt.user( prompt        = "Please give experiment run name - unique, no whitespace",
                                                            default   = default.input$basic$name.exp.run,
                                                            input.con = input.con)

  # Get particle type
  user.input$basic$chosen.particle          <- as.numeric(HX.prompt.user( prompt        = "Choose particles",
                                     choices       = df.particles$nice.particle.name,
                                     default   = default.input$basic$chosen.particle,
                                     input.con = input.con))

  # choose Ripple Filter
  user.input$basic$chosen.rifi.idx             <- as.numeric(HX.prompt.user( prompt        = "Choose ripple filter",
                                        choices       = df.rippleFilter$filter.name,
                                                                  default   = default.input$basic$chosen.rifi.idx,
                                                                  input.con = input.con))
  
  # choose Intensity
  if(df.particles$nice.particle.name[user.input$basic$chosen.particle]%in%c("Protons","Carbon ions","Oxygen ions")){
    user.input$basic$intensity             <- as.numeric(HX.prompt.user( prompt        = "Choose intensity",
                                          choices       = c("System default", paste("I", 3:10, sep = "")),
                                                              default   = default.input$basic$intensity,
                                                              input.con = input.con))
  }else{
    user.input$basic$intensity          <- 0
  }
  
  # Choose spot size
  user.input$basic$spot.distance.mm <- as.numeric(HX.prompt.user( prompt        = "Give fix spot distance in mm, or value <= 0 if distance should be optimized (",
                                                                        default   = default.input$basic$resolution.mm,
                                                                        input.con = input.con))
  return( user.input )
  
  # Choose resolution to optimize within homogenous field
  user.input$basic$resolution.mm          <- as.numeric(HX.prompt.user( prompt        = "Resolution for optimizing homogeneity - in mm, depending e.g. on your detector resolution (if unsure \'1\' is a good guess",
                                                             default   = default.input$basic$resolution.mm,
                                                             input.con = input.con))
  return( user.input )
}