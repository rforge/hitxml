#' Uses weights from SOBP script to create HIT plan
#' Started 2015-08-01, SG

rm(list =ls())
library(HITXML)
library(XML)
library(digest)
library(libamtrack)

#' USER INPUT START

SIS.path           <- "D://04 - Risoe, DKFZ//03 - Methodik//11-20//20 - TRiP//04 - TRiP Basic Data//HIT//03 - TRiP98DATA_HIT-20131120//SIS"
SIS.file           <- "1H_1.1.2009.sis"

SOBP.file          <- "SOBP.dat"
fluence.scaling.factor <- 1.0 # Use if you want to change the fluence in a biologically optimized plan (as changing the
                              # biological dose will not scale the phys. dose / fluence linearily)
particle.name      <- "1H"

name.exp.series    <- "SOBP"
name.exp.run       <- "SOBP_1H"

rifi               <- c("None", "3 mm")[1]
field.shape        <- c("square", "circular")[1]
field.side.size.mm <- 50

focus.no           <- 3

plan.format        <- 2 # plan format always 2 == new

#' USER INPUT END

program.version <- packageDescription("HITXML")$Version
program.date    <- packageDescription("HITXML")$Date

# Translate user input
chosen.particle  <- which(particle.name == c("1H", "4He", "12C", "16O"))
chosen.rifi.idx  <- which(c("None", "3 mm") == rifi)
chosen.field.idx <- which(c("square", "circular") == field.shape)

# Structure to hold all input data
user.input <- list( basic = list( date             = format(Sys.Date(), "%Y%m%d"),
                                  name.exp.series  = name.exp.series,
                                  name.exp.run     = name.exp.series,
                                  chosen.particle  = chosen.particle,
                                  chosen.rifi.idx  = chosen.rifi.idx,        # c("None", "3 mm")
                                  intensity        = 1,
                                  resolution.mm    = 2),
                    
                    fixed = list( time             = format(Sys.time(), "%Y-%m-%dT%H:%M:%S.1111111+02:00"), # Due to POSIX/Windows problem hardcode fractions of second and time zone for timestamp, TODO: Fix by more flexible solution
                                  room.name        = "Room4",
                                  patient.id       = "PT-2004-01",
                                  tumor.name       = "Sacral Chordoma",
                                  patient.birth    = "2004-01-01",
                                  patient.sex      = "M",
                                  therapist.name   = "USER NAME"),
                    
                    IES   = list( ))


  
# Set-up particle data
# TODO: move to external data set, minimum number of particles to input
df.particles <- data.frame( particle.name      = c("PROTON", "ION", "ION", "ION"),
                            nice.particle.name = c("Protons", "Helium-4 ions", "Carbon ions", "Oxygen ions"), 
                            mass               = c(1,4,12,16),
                            charge             = c(1,2,6,8),
                            atomicNumber       = c(1,2,6,8),
                            minParticles       = c(165000, 0, 5000, 0),
                            stringsAsFactors   = FALSE)
df.particles$particle.no <- AT.particle.no.from.Z.and.A(Z = df.particles$charge,
                                                        A = df.particles$mass)$particle.no
df.rippleFilter <- data.frame( filter.name      = c("None", "3mm"),
                               filter.code      = c(254, 3),
                               stringsAsFactors = FALSE)

df.fields       <- data.frame( shape            = c( "square", "circle"),
                               stringsAaFactors = FALSE)

rippleFilter        <- df.rippleFilter$filter.name[user.input$basic$chosen.rifi.idx]
ripplefilter.code   <- df.rippleFilter$filter.code[user.input$basic$chosen.rifi.idx]

# Read SIS data for chosen particle
SIS.data      <- dataSIS(file.name = SIS.file,
                         sis.path  = SIS.path)


# Read SOBP.dat and fill in parameters
SOBP <- read.table( SOBP.file,
                    col.names = c("beam.energy.MeV.u", "V2", "V3", "V4", "fluence.cm2"))[,c(1,5)]
SOBP$beam.energy.MeV.u <- SOBP$beam.energy.MeV.u * 1000

n.IES <- length(SOBP$beam.energy.MeV.u)
for(i in 1:n.IES){
  # i <- 1
  user.input$IES[[i]] <- list( energy.value.MeV.u     = get.clostest.beam.energy.MeV.u(SIS.data, SOBP$beam.energy.MeV.u[i]),
                               chose.idx              = get.IES.index(SIS.data, SOBP$beam.energy.MeV.u[i]),
                               chosen.foc.idx         = focus.no,
                               focus.FWHM.mm          = get.focus.FWHM.mm(SIS.data, SOBP$beam.energy.MeV.u[i], focus.no),
                               field.shape.idx        = chosen.field.idx,
                               fluence.cm2.or.dose.Gy = SOBP$fluence.cm2[i] * fluence.scaling.factor,
                               field.size.mm          = field.side.size.mm,
                               r.min.m                = 15.0)
                         
}

nodes.IES <- list()
df.IES    <- NULL
df.SOBP.feild <- NULL   # for SOBP.dat used in Fluka using the optimized field sized and spots

# Open pdf for report
pdf(paste0(name.exp.run, ".pdf"))
    
for(cur.IES in 1:n.IES){
  # cur.IES <- 1
  # Compute dose from fluence or vice versa
  
  cat("############### Optimizing IES", cur.IES, "#################\n")
  user.input$IES[[cur.IES]]$spill.intensity <- 0
  
  if(user.input$IES[[cur.IES]]$fluence.cm2.or.dose.Gy < 0.0){
    dose.Gy     <- -1.0 * user.input$IES[[cur.IES]]$fluence.cm2.or.dose.Gy
    fluence.cm2 <- AT.fluence.cm2.from.dose.Gy( E.MeV.u      = user.input$IES[[cur.IES]]$energy.value.MeV.u,
                                                D.Gy         = dose.Gy,
                                                particle.no  = df.particles$particle.no[user.input$basic$chosen.particle],
                                                material.no  = 1,                          # i.e. liquid water
                                                stopping.power.source.no = 3)$fluence.cm2  # use ICRU49/73 data
  }else{                                             
    fluence.cm2 <- user.input$IES[[cur.IES]]$fluence.cm2.or.dose.Gy
    dose.Gy     <- AT.dose.Gy.from.fluence.cm2( E.MeV.u      = user.input$IES[[cur.IES]]$energy.value.MeV.u,
                                                fluence.cm2  = fluence.cm2,
                                                particle.no  = df.particles$particle.no[user.input$basic$chosen.particle],
                                                material.no  = 1,                      # i.e. liquid water
                                                stopping.power.source.no = 3)$dose.Gy  # use ICRU49/73 data
  }   
  
  # Prepare start and field parameters
  if(df.fields$shape[user.input$IES[[cur.IES]]$field.shape.idx] == "circle"){
    par.start    <- c(user.input$IES[[cur.IES]]$focus.FWHM.mm / 2, 
                      user.input$IES[[cur.IES]]$focus.FWHM.mm / 3)
    field.par    <- c(user.input$IES[[cur.IES]]$field.size.mm, 
                      user.input$IES[[cur.IES]]$r.min.mm)
  }
  if(df.fields$shape[user.input$IES[[cur.IES]]$field.shape.idx] == "square"){
    par.start    <- user.input$IES[[cur.IES]]$focus.FWHM.mm * 0.99
    field.par    <- user.input$IES[[cur.IES]]$field.size.mm
  }
  
  # Optimize field
  optim.beam.spot.grid <- HX.optimize.field( field.shape   = df.fields$shape[user.input$IES[[cur.IES]]$field.shape.idx],
                                             par.start     = par.start,
                                             focus.FWHM.mm = user.input$IES[[cur.IES]]$focus.FWHM.mm,
                                             field.par     = field.par, 
                                             fluence.cm2   = fluence.cm2,
                                             N.min         = df.particles$minParticles[user.input$basic$chosen.particle],
                                             resolution.mm = user.input$basic$resolution.mm,
                                             n.IES         = cur.IES,
                                             plot          = TRUE)$beam.spot.grid
  
  ## Translate beam spot grid into XLM structure
  voxel.node.list <- list()
  for (i in 1:nrow(optim.beam.spot.grid)){
    # i <- 1
    voxel.node.list[[i]] <- xmlNode("Voxel", 
                                    attrs = c( x         = optim.beam.spot.grid$x.mm[i], 
                                               y         = optim.beam.spot.grid$y.mm[i], 
                                               particles = optim.beam.spot.grid$N.particles[i]))
  }
  
  nodes.IES[[cur.IES]] <- xmlNode("IES", 
                                attrs        = c(number          = as.numeric(cur.IES),
                                                 energy          = as.numeric(user.input$IES[[cur.IES]]$energy.value.MeV.u),
                                                 focus           = as.numeric(user.input$IES[[cur.IES]]$focus.FWHM.mm)),
                                .children    = voxel.node.list)
  
  # Add to IES overview data frame
  df.IES <- rbind(df.IES,
                  data.frame( IES.no        = cur.IES,
                              E.MeV.u       = user.input$IES[[cur.IES]]$energy.value.MeV.u,
                              focus.FWHM.mm = user.input$IES[[cur.IES]]$focus.FWHM.mm,
                              fluence.cm2   = fluence.cm2,
                              dose.Gy       = dose.Gy,
                              n.spots       = nrow(optim.beam.spot.grid), # new field added to table (Yukihara, Oct. 2014)
                              part.spot     = optim.beam.spot.grid$N.particles[1] # new field added to table (Yukihara, Oct. 2014)
                  ))
  
  
  df.SOBP.feild <- rbind(df.SOBP.feild,
                       
                       data.frame(
                         E.GeV = rep(user.input$IES[[cur.IES]]$energy.value.MeV.u/1000,nrow(optim.beam.spot.grid)),
                         x.cm=optim.beam.spot.grid$x.mm/10,
                         y.cm=optim.beam.spot.grid$y.mm/10,
                         FWHM.cm=rep(user.input$IES[[cur.IES]]$focus.FWHM.mm/10,nrow(optim.beam.spot.grid)),
                         N.particles =optim.beam.spot.grid$N.particles
                       )
                       
                       
  )
}

row.names(df.IES) <- 1:nrow(df.IES)
clan.txtplot(df.IES)

# Write plan documentation in detail (Yukihara, Oct. 2014)
documentation <- paste("\n\nHIT_XML detailed plan documentation\n\n",
                       "Plan name:",
                       user.input$basic$date, "/", 
                       user.input$basic$name.exp.series, "/", 
                       user.input$basic$name.exp.run, "\n",
                       "Projectile: ",
                       df.particles$nice.particle.name[user.input$basic$chosen.particle],
                       "\n",
                       sep = "")
for (i in 1:length(user.input$basic)){
  documentation <- paste(documentation,
                         names(user.input$basic)[i]," \t\t\t = \t\t\t",as.character(user.input$basic[i]),"\n",
                         sep="")}

for (i in 1:length(user.input$IES)){
  documentation <- paste(documentation, "\nIES no.", i, "\n", sep = "")
  for (j in 1:length(user.input$IES[[i]])){
    documentation <- paste(documentation,
                           names(user.input$IES[[i]][j])," = ",as.character(user.input$IES[[i]][[j]]),"\n",
                           sep = "")
  }
  
}


clan.txtplot(documentation)
dev.off()



first.name.exp.run <- name.exp.run

# Write request
file.name.plan    <- paste(user.input$basic$name.exp.run, ".xml", sep = "")
file.name.request <- paste("req_", file.name.plan, sep = "")
file.name.result  <- paste(user.input$basic$name.exp.run, "_result.xml", sep = "")
file.name.PBR     <- paste(user.input$basic$name.exp.run, "_PBR.xml", sep = "")
file.name.MBR     <- paste(user.input$basic$name.exp.run, "_MBR.xml", sep = "")

path.main   <- "."
path.plans  <- paste(path.main, user.input$basic$date, "\\", user.input$basic$name.exp.series, "\\", sep = "")
path.result <- paste(path.plans, "OUTPUT\\", sep = "")
path.PBR    <- paste(path.result, "PhysicalBeamRecordFile\\", file.name.PBR, sep = "")
path.MBR    <- paste(path.result, "MachineBeamRecordFile\\", file.name.MBR, sep = "")
path.save   <- path.main

program.context <- paste(path.plans, file.name.plan, sep = "")
md5.checksum    <- digest(program.context, algo = "md5")

node.GroupParameter <- xmlNode(  "GroupParameter", attrs = c(name = "ExecuteSFP"),
                               .children = list(xmlNode("Parameter", 
                                                        attrs = c(name = "ProgramID"), 
                                                        value = "TCsPerformIrradiation"),
                                                xmlNode("Parameter", 
                                                        attrs = c(name = "ProgramContext"), 
                                                        value = program.context),
                                                xmlNode("Parameter", 
                                                        attrs = c(name = "ProgramContextMd5CheckSum"), 
                                                        value = md5.checksum),
                                                xmlNode("Parameter", 
                                                        attrs = c(name = "ProgramResult"), 
                                                        value = (paste(path.result, file.name.result, sep = ""))),
                                                xmlNode("Parameter", 
                                                        attrs = c(name = "Timeout"), 
                                                        value = "3000000")))

node.Message        <- xmlNode( "Message", 
                              attrs = c( type      = "SET", 
                                         timestamp = user.input$basic$time, 
                                         device    = "SFPGM"), 
                              .children = list(node.GroupParameter))

## Write plan

## The "Housekeeping", "Specifications", and "UI" nodes are not
## necessary and not used. Their generation code is kept
## in the svn of <0.6.2 versions. Maybe it will be useful one day.

node.PBR <- xmlNode("Parameter", 
                  attrs     = c(Name = "PhysicalBeamRecordFile"), 
                  value     = path.PBR)
node.MBR <- xmlNode("Parameter", 
                  attrs     = c(Name = "MachineBeamRecordFile"), 
                  value     = path.MBR)

node.Spill <- xmlNode("Parameter", 
                    attrs = c(Name = "FullSpill"), 
                    value = "false")

node.SpillIntensity <- xmlNode("Parameter", 
                             attrs = c(Name = "SpillIntensity"), 
                             value = sprintf("%2.1e", user.input$IES[[1]]$spill.intensity))

node.PosCorrLoop <- xmlNode("Parameter", 
                          attrs = c(Name = "EnablePosCorrectionLoop"), 
                          value = "true")

node.DosemeterName <- xmlNode("Parameter", 
                            attrs = c(Name = "DosemeterDeviceName"), 
                            value = "none")

node.DosemeterCommand <- xmlNode("Parameter", 
                               attrs = c(Name = "DosemeterCommand"), 
                               value = "none")

node.Channel <- xmlNode("Parameter", 
                      attrs = c(Name = "Channel"), 
                      value = "All")

node.EnableIC1 <- xmlNode("Parameter", 
                        attrs = c(Name = "EnableIC1"), 
                        value = "true")

node.EnableIC2 <- xmlNode("Parameter", 
                        attrs = c(Name = "EnableIC2"), 
                        value = "true")

node.EnableIM <- xmlNode("Parameter", 
                       attrs = c(Name = "EnableIM"), 
                       value = "true")

node.RstFormat <- xmlNode( "RstFormat", 
                         value = "PT_2004")

node.Patient <-  xmlNode("Patient", 
                       attrs = c( id        = user.input$fixed$patient.id, 
                                  name      = user.input$fixed$tumor.name, 
                                  sex       = user.input$fixed$patient.sex, 
                                  birthDate = user.input$fixed$patient.birth))

node.TxInitiation <-  xmlNode("TxInitiation", 
                            attrs = c( therapist = user.input$fixed$therapist.name, 
                                       dateTime  = "2007-01-23T13:52:27.2343750+01:00"))

node.BAMS <- xmlNode( "BAMS", 
                    attrs = c( rippleFilter         = ripplefilter.code , 
                               rangeShifter         = "3", 
                               rangeShifterDistance = "20"))
node.table <- xmlNode( "TxTable", 
                     attrs = c( roll            = "0", 
                                pitch           = "0", 
                                lateral         = "150", 
                                longitudinal    = "-400", 
                                isocentricAngle = "0", 
                                vertical        = "-100.5"))
node.gantry <- xmlNode( "Gantry", 
                      attrs = c( angle = "90"))

tmp.projectile   <- df.particles$particle.name[user.input$basic$chosen.particle]
tmp.charge       <- df.particles$charge[user.input$basic$chosen.particle]
tmp.mass         <- df.particles$mass[user.input$basic$chosen.particle]
tmp.atomicNumber <- df.particles$atomicNumber[user.input$basic$chosen.particle]

if(tmp.projectile == "PROTON"){    # Funny enough, but as from Jun06 mass/charge/atomic no for protons are "" rather than 1
tmp.charge		<- ""
tmp.mass		<- ""
tmp.atomicNumber<- ""
}


node.TxRoom <- xmlNode("TxRoom", 
                     attrs = c( name         = user.input$fixed$room.name, 
                                projectile   = tmp.projectile, 
                                charge       = tmp.charge, 
                                mass         = tmp.mass, 
                                atomicNumber = tmp.atomicNumber))

# REMARK: The beam UID could be changed in principle but it
# likely to cause problems as the same beam UID must be provided manually
# to cancel beam requests in the queueing system, so just leave it
# as it is
node.Beam <- xmlNode("Beam", 
                   attrs     = c(uid = "bee035c5-03f6-4e9c-94b9-31f0fc484db1"),                           
                   .children = c(list( node.RstFormat,
                                       node.Patient,
                                       node.TxInitiation,
                                       node.TxRoom,
                                       node.BAMS,
                                       node.table,
                                       node.gantry),
                                 nodes.IES))

node.PTTxPlan <- xmlNode( "PTTxPlan", 
                        .children = list(node.Beam))

node.BeamPlan <- xmlNode( "Parameter", 
                        attrs = c(Name = "BeamPlan"),
                        .children = list(node.PTTxPlan))

if(user.input$IES[[1]]$spill.intensity == 0){
node.SfpParameters <- xmlNode( "SfpParameters", 
                               .children = list( node.PBR, 
                                                 node.MBR, 
                                                 node.Spill, 
                                                 node.PosCorrLoop, 
                                                 node.DosemeterName, 
                                                 node.DosemeterCommand, 
                                                 node.Channel,
                                                 node.EnableIC1,
                                                 node.EnableIC2, 
                                                 node.EnableIM,
                                                 node.BeamPlan))
}else{
node.SfpParameters <- xmlNode( "SfpParameters", 
                               .children = list( node.PBR, 
                                                 node.MBR, 
                                                 node.Spill, 
                                                 node.SpillIntensity, 
                                                 node.PosCorrLoop, 
                                                 node.DosemeterName, 
                                                 node.DosemeterCommand, 
                                                 node.Channel,
                                                 node.EnableIC1,
                                                 node.EnableIC2, 
                                                 node.EnableIM,
                                                 node.BeamPlan))
}



# It seems we do not need the first three node for a working plan
node.WorkflowContext <- xmlNode( "WorkflowContext",
                               .children = list( #node.Housekeeping,
                                 #node.Specifications,
                                 #node.UIs,
                                 node.SfpParameters))

# 	# THIS IS THE NEW PLAN FORMAT
# 	node.PTTxPlanMd5 <- xmlNode("PTTxPlanMd5",
# 								.children = node.PTTxPlan,
# 								attrs = c( md5         						= "noMD5")) 
# 									       xmlns:xsi   						= "http://www.w3.org/2001/XMLSchema-instance", 
# 									       xsi:noNamespaceSchemaLocation    = "RTT-PT-Plan.xsd"))

# THIS IS THE NEW PLAN FORMAT - with small fixes (Yukihara, Oct. 2014)
node.PTTxPlanMd5 <- xmlNode("PTTxPlanMd5",
                          node.PTTxPlan,
                          attrs = c( md5 = "noMD5",
                                     "xmlns:xsi" = "http://www.w3.org/2001/XMLSchema-instance", 
                                     "xsi:noNamespaceSchemaLocation" = "RTT-PT-Plan.xsd"))

# SAVE
dir.create(path.save,
         recursive = TRUE,
         showWarnings = FALSE)


prefix <- paste("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
              "<!-- IRRADIATION PLAN CREATED WITH HIT_XML ",
              program.version,
              " (",
              program.date,
              ") -->\n",
              sep = "")

if(plan.format == 1){
saveXML( node.Message, 
         file = file.path(paste(path.save, "/", file.name.request, sep = "")),
         prefix = prefix)

saveXML( node.WorkflowContext, 
         file = file.path(paste(path.save, "/", file.name.plan, sep = "")),
         prefix = prefix)

cat("\nSaved request file ",
    file.name.request,
    " and plan file ",
    file.name.plan,
    " to location ",
    path.save,
    "\n", sep = "")
}else{

# Save plan in the new format, according to Steffen's format (Yukihara, Oct. 2014)
saveXML( node.PTTxPlanMd5, 
         file = file.path(paste(path.save, "/", file.name.plan, sep = "")),
         prefix = prefix)

cat("\nSaved plan file ",
    file.name.plan,
    " to location ",
    path.save,
    "\n",
    sep = "")
}




# saving the SOBP.dat version of the xml plan for FLUKA simulations

write.table(df.SOBP.feild,
            file = "SOBP-optimized.dat",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol       = "\r\n" )





