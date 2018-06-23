#!/usr/bin/Rscript

############################################
# Script to create HIT irradiation plans for
# simple, homogenous fields
#
# E0408, DKFZ Heidelberg
# Oct2015
# Maintainer: s.greilich@dkfz.de
############################################
# Please see README for more information
############################################

############################################
# TODOs
############################################
# Fine tuning parameters
# Time estimate (completed in Oct. 2014 - Yukihara)
############################################

# Modifications implemented (Oct. 2014, Yukihara) 
# - small fixes to the code to allow it to generate multiple IES with
#   different energies. 
# - Modified to allow the repetition of IES without repeating the optimization. 
# - The code now saves by default the plans in the new format (tested successfully
#   at HIT)
#
# Known problems
# - The function that gathers user input does not contain He ion and, therefore,
#   does not ask for intensity. As a result, time estimate does not work. The plan
#   should work without problems, since intensity information is not included
#   in the plans.


# Clear workspace, set version
rm(list =ls())

# load package
library(HITXML)

program.version <- packageDescription("HITXML")$Version
program.date    <- packageDescription("HITXML")$Date

# Structure to hold all input data
# This approach allows to run the script interactively, too,
# e.g. for debugging. Normal mode is non-interactively
# by Rscript from command line.
user.input <- list( basic = list( date             = 20160920,
                                  name.exp.series  = "anTest",
                                  name.exp.run     = "ab0001Test",
                                  chosen.particle  = 1,        # c("Protons", "Helium-4 ions", "Carbon ions", "Oxygen ions")
                                  chosen.rifi.idx  = 1,        # c("None", "3 mm")
                                  intensity        = 1,
                                  resolution.mm    = 1,
                                  spot.distance.mm = -1),      # If <= 0, do spot distance optimization, otherwise use given value
                    
                    fixed = list( time             = format(Sys.time(), "%Y-%m-%dT%H:%M:%S.1111111+02:00"), # Due to POSIX/Windows problem hardcode fractions of second and time zone for timestamp, TODO: Fix by more flexible solution
                                  room.name        = "Room4",
                                  patient.id       = "PT-2004-01",
                                  tumor.name       = "Sacral Chordoma",
                                  patient.birth    = "2004-01-01",
                                  patient.sex      = "M",
                                  therapist.name   = "USER NAME"),
                    
                    IES   = list( list( energy.value.MeV.u     = 50.6,       # value will be adjusted to closest available in libC
                                        chosen.idx             = NA,         # if given, overrides energy value
                                        focus.FWHM.mm          = 0.0,        # value will be adjusted to closest available in libC
                                        chosen.foc.idx         = 3,          # if given, overrides focus value
                                        field.shape.idx        = 1,          # c("square", "circular", "shells")
                                        fluence.cm2.or.dose.Gy = 2e5,
                                        field.size.mm          = 50,
                                        r.min.mm               = 0.0)))

# Script run interactively or non-interactively, set connection accordingly
interactive.run.mode <- interactive()
if(!interactive.run.mode){
  input.con <-file("stdin")
}else{
  input.con <-stdin()
}

# Read defaults
default.input <- HX.read.defaults()

# Say hello to user
welcome.message <- paste( "\n",
	 "#####################################################################\n",
	 "#####################################################################\n",
	 "## THIS IS HIT_XML\n## E0408, DKFZ Heidelberg\n## Version ", 
     program.version, 
	 " as of ", 
	 program.date,
	 "\n",
	 "#####################################################################\n",
     sep = "")

message(welcome.message)

message("#####################################################################\n\n")

#################################
## Set up data frames for choices
df.particles <- data.frame( particle.name      = c("PROTON", "ION", "ION", "ION"),
                            nice.particle.name = c("Protons", "Helium ions", "Carbon ions", "Oxygen ions"), 
                            mass               = c(1,4,12,16),
                            charge             = c(1,2,6,8),
                            atomicNumber       = c(1,2,6,8),
                            minParticles       = c(165000, 0, 5000, 0),
                            libC.file.name     = c( "LIBC HIT Version 1.4 p V11.txt", 
                                                    "LIBC HIT Version 1.5 He4 V1.txt",
                                                    "LIBC HIT Version 1.4 C12 V8.txt",
                                                    "LIBC HIT Version 1.5 O16 V2.txt"),
                            stringsAsFactors   = FALSE)

df.particles$particle.no <- libamtrack::AT.particle.no.from.Z.and.A( Z = df.particles$charge,
                                                                     A = df.particles$mass)$particle.no

df.rippleFilter <- data.frame( filter.name      = c("None", "3mm"),
                               filter.code      = c(254, 3),
                               stringsAsFactors = FALSE)

df.fields       <- data.frame( shape            = c( "square", "circle", "shells"),
                               stringsAaFactors = FALSE)

# TODO: Add variables for field shape here

#######################################
## Get user inputs (if non-interactive)

# There is currently a problem wiht the function HX.get.user.input.basic.
# The function does not include He ions in a list of possible particles,
# so it is necessary to update it.
# (Yukihara, 30 Oct 14)

if(!interactive.run.mode){
  user.input <- HX.get.user.input.basic(user.input      = user.input,
                                        default.input   = default.input,
                                        df.particles    = df.particles, 
                                        df.rippleFilter = df.rippleFilter,
                                        input.con       = input.con)  
} 

# Translate ripple filter code
rippleFilter        <- df.rippleFilter$filter.name[user.input$basic$chosen.rifi.idx]
ripplefilter.code   <- df.rippleFilter$filter.code[user.input$basic$chosen.rifi.idx]

# load libC for chosen particle
df.libc.HIT   <- read.table( system.file( "extdata", "libC",
                                          df.particles$libC.file.name[user.input$basic$chosen.particle],
                                          package = "HITXML"), 
							 sep        = "", 
							 col.names = c( "energy.no", "energy.MeV.u", 
					   						paste("F", 1:6,  ".mm",  sep = ""),
					                        paste("I", 1:15, ".1.s", sep = "")))

# Simple path (Yukihara, Oct. 2014)
path.save   <- file.path( gsub( "\\", 
                                "/", 
                                paste( "./", 
                                       user.input$basic$date,
                                       sep = ""),
                                fixed = TRUE))
dir.create(path.save,
           recursive = TRUE,
           showWarnings = FALSE)

if(!interactive.run.mode){
  pdf(file = paste(path.save,"/",user.input$basic$name.exp.run, ".pdf", sep = ""))
}

title.doc <- paste(welcome.message,
				   "\n\nHIT_XML plan documentation\n\n",
				   "Plan name:",
           user.input$basic$date, "/", 
           user.input$basic$name.exp.series, "/", 
           user.input$basic$name.exp.run, "\n",
				   "Projectile: ",
				   df.particles$nice.particle.name[user.input$basic$chosen.particle],
				   "\n",
				   sep = "")

clan.txtplot(title.doc)

n.IES     <- 1
nodes.IES <- list()
df.IES    <- NULL

repeat{
	cat("\n#####################################################################\n",
		"Generating iso-energy slice (IES) no. ", n.IES, ".\n", 
		"#####################################################################\n\n",
		sep = "")
	
	if(!interactive.run.mode){
	  user.input <- HX.get.user.input.IES(  n.IES           = n.IES,
                                          user.input      = user.input,
	                                        default.input   = default.input,
	                                        df.particles    = df.particles,
                                          df.fields       = df.fields,
	                                        df.libc.HIT     = df.libc.HIT,
                                          input.con       = input.con)
	}else{
    user.input <- HX.correct.user.input(user.input     = user.input, 
                                        n.IES          = n.IES, 
                                        df.libc.HIT    = df.libc.HIT)
	}

  cat("Compute dose from fluence or vice versa\n")
	# Compute dose from fluence or vice versa
  if(user.input$IES[[n.IES]]$fluence.cm2.or.dose.Gy < 0.0){
	  dose.Gy     <- -1.0 * user.input$IES[[n.IES]]$fluence.cm2.or.dose.Gy
	  fluence.cm2 <- libamtrack::AT.fluence.cm2.from.dose.Gy( E.MeV.u      = user.input$IES[[n.IES]]$energy.value.MeV.u,
            	                                              D.Gy         = dose.Gy,
            	                                              particle.no  = df.particles$particle.no[user.input$basic$chosen.particle],
            	                                              material.no  = 1,                          # i.e. liquid water
            	                                              stopping.power.source.no = 3)$fluence.cm2  # use ICRU49/73 data
	}else{                                             
	  fluence.cm2 <- user.input$IES[[n.IES]]$fluence.cm2.or.dose.Gy
	  dose.Gy     <- libamtrack::AT.dose.Gy.from.fluence.cm2( E.MeV.u      = user.input$IES[[n.IES]]$energy.value.MeV.u,
            	                                              fluence.cm2  = fluence.cm2,
            	                                              particle.no  = df.particles$particle.no[user.input$basic$chosen.particle],
            	                                              material.no  = 1,                      # i.e. liquid water
            	                                              stopping.power.source.no = 3)$dose.Gy  # use ICRU49/73 data
	}   
	
  cat("Prepare start and field parameters\n")
  # Prepare start and field parameters
  if(df.fields$shape[user.input$IES[[n.IES]]$field.shape.idx] == "shells"){
    par.start    <- user.input$IES[[n.IES]]$focus.FWHM.mm * 0.99
    field.par    <- user.input$IES[[n.IES]]$field.size.mm
  }
  if(df.fields$shape[user.input$IES[[n.IES]]$field.shape.idx] == "circle"){
	   par.start    <- c(user.input$IES[[n.IES]]$focus.FWHM.mm / 2, 
                       user.input$IES[[n.IES]]$focus.FWHM.mm / 3)
		 field.par    <- c(user.input$IES[[n.IES]]$field.size.mm, 
		                   user.input$IES[[n.IES]]$r.min.mm)
	}
	if(df.fields$shape[user.input$IES[[n.IES]]$field.shape.idx] == "square"){
		 par.start    <- user.input$IES[[n.IES]]$focus.FWHM.mm * 0.99
		 field.par    <- user.input$IES[[n.IES]]$field.size.mm
	}
  
  cat("Optimize field\n")
  # Optimize field
	optim.beam.spot.grid <- HX.optimize.field( field.shape   = df.fields$shape[user.input$IES[[n.IES]]$field.shape.idx],
		                                         par.start     = par.start,
												                     focus.FWHM.mm = user.input$IES[[n.IES]]$focus.FWHM.mm,
												                     field.par     = field.par, 
												                     fluence.cm2   = fluence.cm2,
												                     N.min         = df.particles$minParticles[user.input$basic$chosen.particle],
												                     resolution.mm = user.input$basic$resolution.mm,
												                     spot.distance.mm = user.input$basic$spot.distance.mm,
												                     n.IES         = n.IES,
												                     plot          = TRUE)$beam.spot.grid

	# TODO: Add time estimate for plan
	# TODO: Add robustness analysis by MC (see older versions in svn for autofocus routine)
											   
	## Translate beam spot grid into XLM structure
	voxel.node.list <- list()
	for (i in 1:nrow(optim.beam.spot.grid)){
		 # i <- 1
		 voxel.node.list[[i]] <- XML::xmlNode("Voxel", 
                    										   attrs = c( x         = optim.beam.spot.grid$x.mm[i], 
                  	  												        y         = optim.beam.spot.grid$y.mm[i], 
                  		  											        particles = optim.beam.spot.grid$N.particles[i]))
	}

	nodes.IES[[n.IES]] <- XML::xmlNode("IES", 
                					    	      attrs        = c(number          = as.numeric(n.IES),
                    								   	               energy          = as.numeric(user.input$IES[[n.IES]]$energy.value.MeV.u),
                		           									       focus           = as.numeric(user.input$IES[[n.IES]]$focus.FWHM.mm)),
                						                           .children    = voxel.node.list)

	# Add to IES overview data frame
	df.IES <- rbind(df.IES,
					data.frame( IES.no        = n.IES,
								E.MeV.u       = user.input$IES[[n.IES]]$energy.value.MeV.u,
								focus.FWHM.mm = user.input$IES[[n.IES]]$focus.FWHM.mm,
								fluence.cm2   = fluence.cm2,
								dose.Gy       = dose.Gy,
                n.spots       = nrow(optim.beam.spot.grid), # new field added to table (Yukihara, Oct. 2014)
                part.spot     = optim.beam.spot.grid$N.particles[1] # new field added to table (Yukihara, Oct. 2014)
                  ))
  
  #Repeat IES? (Yukihara, Oct 2014)
  if(!interactive.run.mode){
    cat("\nRepeat IES? Enter TOTAL number of identical IES [press return for 1]?\n-->",sep="")
    answer<-readLines(con=input.con,n=1)
    n.repeats=as.numeric(answer)-1
    if(!is.na(n.repeats) & n.repeats>0){
      cat("\nCreating ",n.repeats," IES's, for a total of ",n.repeats+1," identical IES's.\n",sep="")
      for(i in (n.IES+1):(n.IES+n.repeats)){
        #repeat user input  
        user.input$IES[[i]]=user.input$IES[[i-1]]
      
        # Repeat node
        nodes.IES[[i]] <- XML::xmlNode("IES", 
                                       attrs        = c(number          = as.numeric(i),
                                                        energy          = as.numeric(user.input$IES[[i]]$energy.value.MeV.u),
                                                        focus           = as.numeric(user.input$IES[[i]]$focus.FWHM.mm)),
                                       .children    = voxel.node.list)
      
        # Add to IES overview data frame
        df.IES <- rbind(df.IES,
                        data.frame( IES.no        = i,
                                    E.MeV.u       = user.input$IES[[i]]$energy.value.MeV.u,
                                    focus.FWHM.mm = user.input$IES[[i]]$focus.FWHM.mm,
                                    fluence.cm2   = fluence.cm2,
                                    dose.Gy       = dose.Gy,
                                    n.spots       = nrow(optim.beam.spot.grid), # new field added to table (Yukihara, Oct. 2014)
                                    part.spot     = optim.beam.spot.grid$N.particles[1] # new field added to table (Yukihara, Oct. 2014)
                        ))
      }
      n.IES=n.IES+n.repeats
    }
  }
	
  # Another IES?
	if(!interactive.run.mode){
    cat( "\nGenerate another IES [y/n (default, press return)]?\n--> ",
			   sep = "")
	  answer <- readLines(con     = input.con,
	                      n       = 1)
	  if(answer != 'y'){
		  break
	  }
    user.input$IES=c(user.input$IES,list(list())) # Increase dimention to add IES (Yukihara, Oct. 2014)
  }else{
    break
  }

	n.IES <- n.IES + 1
} # repeat IES

row.names(df.IES) <- 1:nrow(df.IES)
clan.txtplot(df.IES)

# Estimate irradiation time (Yukihara, Oct. 2014)
total.particles=sum(df.IES$n.spots*df.IES$part.spot)
duty.cycle=0.4

if(user.input$basic$intensity == 0){
  spill.intensity <- df.libc.HIT[1, 8+10]
  cat("\n\n=====================================================")
  cat("\nIntensity not defined.")
  cat("\n=====================================================\n")
}else{
  spill.intensity <- df.libc.HIT[1,8 + user.input$basic$intensity]
}
cat("\nCalculating irradiation time with intensity = ",spill.intensity," particles/s")

# spill.intensity=user.input$IES[[1]]$spill.intensity
total.time.s=total.particles/spill.intensity/duty.cycle

time.estimation.message=paste(
  "\nTime of irradiation (estimate) and total dose",
  "\n=================================================================",
  "\n Total number of particles     =",total.particles,
  "\n Spill intensity (particles/s) =",spill.intensity,
  "\n Duty cycle                    =",duty.cycle,
  "\n=================================================================",
  "\n Estimated time (s)            =",total.time.s,
  "\n Estimate time (min),          =",total.time.s/60,
  "\n=================================================================",
  "\n Total dose (Gy)            =",sum(df.IES$dose.Gy),
  "\n=================================================================\n"
  )
cat(time.estimation.message)
clan.txtplot(time.estimation.message)

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


if(!interactive.run.mode){
  dev.off()
}

# new or old plan format? Non-interactive --> always 2 == new
if(!interactive.run.mode){
    plan.format <- as.numeric( HX.prompt.user( prompt        = "\nWrite plan in new or old format?",
                                               choices       = c("Old", "New"),
                                               default       = 2,
                                               input.con     = input.con))
}else{
	plan.format <- 2
}

first.name.exp.run <- user.input$basic$name.exp.run

repeat{
	# Write back defaults
	HX.save.defaults(user.input)

	# Write request
	file.name.plan    <- paste(user.input$basic$name.exp.run, ".xml", sep = "")
	file.name.request <- paste("req_", file.name.plan, sep = "")
	file.name.result  <- paste(user.input$basic$name.exp.run, "_result.xml", sep = "")
	file.name.PBR     <- paste(user.input$basic$name.exp.run, "_PBR.xml", sep = "")
	file.name.MBR     <- paste(user.input$basic$name.exp.run, "_MBR.xml", sep = "")

	path.main   <- "E:\\Rtt-Pt-Sim-WorkDir\\HIT-Exp\\"
	path.plans  <- paste(path.main, user.input$basic$date, "\\", user.input$basic$name.exp.series, "\\", sep = "")
	path.result <- paste(path.plans, "OUTPUT\\", sep = "")
	path.PBR    <- paste(path.result, "PhysicalBeamRecordFile\\", file.name.PBR, sep = "")
	path.MBR    <- paste(path.result, "MachineBeamRecordFile\\", file.name.MBR, sep = "")

#   # Old directory structure from Steffen (removed by Yukihara, Oct. 2014)
# 	path.save   <- file.path( gsub( "\\", 
# 									"/", 
# 										paste( "./Rtt-Pt-Sim-WorkDir/HIT-Exp/", 
# 										       user.input$basic$date, 
# 											   "/", 
# 										       user.input$basic$name.exp.series, 
# 											   sep = ""),
# 									fixed = TRUE))
  


	program.context <- paste(path.plans, file.name.plan, sep = "")
	md5.checksum    <- digest::digest(program.context, algo = "md5")

	node.GroupParameter <- XML::xmlNode(  "GroupParameter", 
	                                      attrs = c(name = "ExecuteSFP"),
										                    .children = list(XML::xmlNode("Parameter", 
																        attrs = c(name = "ProgramID"), 
																        value = "TCsPerformIrradiation"),
												 XML::xmlNode("Parameter", 
      																 attrs = c(name = "ProgramContext"), 
      																 value = program.context),
												 XML::xmlNode("Parameter", 
      																 attrs = c(name = "ProgramContextMd5CheckSum"), 
      																 value = md5.checksum),
												 XML::xmlNode("Parameter", 
      																 attrs = c(name = "ProgramResult"), 
      																 value = (paste(path.result, file.name.result, sep = ""))),
												 XML::xmlNode("Parameter", 
            													 attrs = c(name = "Timeout"), 
            													 value = "3000000")))

	node.Message        <- XML::xmlNode( "Message", 
									                     attrs = c( type      = "SET", 
                          											  timestamp = user.input$basic$time, 
                          											  device    = "SFPGM"), 
									                     .children = list(node.GroupParameter))

	## Write plan

	## The "Housekeeping", "Specifications", and "UI" nodes are not
	## necessary and not used. Their generation code is kept
	## in the svn of <0.6.2 versions. Maybe it will be useful one day.

	node.PBR <- XML::xmlNode( "Parameter", 
                						attrs     = c(Name = "PhysicalBeamRecordFile"), 
                						value     = path.PBR)
	node.MBR <- XML::xmlNode( "Parameter", 
                						attrs     = c(Name = "MachineBeamRecordFile"), 
                						value     = path.MBR)

	node.Spill <- XML::xmlNode( "Parameter", 
                						  attrs = c(Name = "FullSpill"), 
                						  value = "false")

	node.SpillIntensity <- XML::xmlNode("Parameter", 
                        						  attrs = c(Name = "SpillIntensity"), 
                        						  value = sprintf("%2.1e", user.input$IES[[1]]$spill.intensity))

	node.PosCorrLoop <- XML::xmlNode( "Parameter", 
                    								attrs = c(Name = "EnablePosCorrectionLoop"), 
                    								value = "true")

	node.DosemeterName <- XML::xmlNode("Parameter", 
                  									 attrs = c(Name = "DosemeterDeviceName"), 
                  									 value = "none")

	node.DosemeterCommand <- XML::xmlNode( "Parameter", 
                      									 attrs = c(Name = "DosemeterCommand"), 
                      									 value = "none")

	node.Channel <- XML::xmlNode( "Parameter", 
                  							attrs = c(Name = "Channel"), 
                  							value = "All")

	node.EnableIC1 <- XML::xmlNode( "Parameter", 
                  							  attrs = c(Name = "EnableIC1"), 
                  							  value = "true")

	node.EnableIC2 <- XML::xmlNode( "Parameter", 
                  							  attrs = c(Name = "EnableIC2"), 
                  							  value = "true")

	node.EnableIM <- XML::xmlNode( "Parameter", 
                  							 attrs = c(Name = "EnableIM"), 
                  							 value = "true")

	node.RstFormat <- XML::xmlNode( "RstFormat", 
                							   value = "PT_2004")

	node.Patient <-  XML::xmlNode("Patient", 
                							 attrs = c( id        = user.input$fixed$patient.id, 
                										name      = user.input$fixed$tumor.name, 
                										sex       = user.input$fixed$patient.sex, 
                										birthDate = user.input$fixed$patient.birth))

	node.TxInitiation <-  XML::xmlNode("TxInitiation", 
								  attrs = c( therapist = user.input$fixed$therapist.name, 
											 dateTime  = "2007-01-23T13:52:27.2343750+01:00"))

	node.BAMS <- XML::xmlNode( "BAMS", 
						  attrs = c( rippleFilter         = ripplefilter.code , 
									 rangeShifter         = "3", 
									 rangeShifterDistance = "20"))
	node.table <- XML::xmlNode( "TxTable", 
						   attrs = c( roll            = "0", 
									  pitch           = "0", 
									  lateral         = "150", 
									  longitudinal    = "-400", 
									  isocentricAngle = "0", 
									  vertical        = "-100.5"))
	node.gantry <- XML::xmlNode( "Gantry", 
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
	
	
	node.TxRoom <- XML::xmlNode("TxRoom", 
						   attrs = c( name         = user.input$fixed$room.name, 
									  projectile   = tmp.projectile, 
									  charge       = tmp.charge, 
									  mass         = tmp.mass, 
									  atomicNumber = tmp.atomicNumber))

	# REMARK: The beam UID could be changed in principle but it
	# likely to cause problems as the same beam UID must be provided manually
	# to cancel beam requests in the queueing system, so just leave it
	# as it is
	node.Beam <- XML::xmlNode("Beam", 
						  attrs     = c(uid = "bee035c5-03f6-4e9c-94b9-31f0fc484db1"),                           
						  .children = c(list( node.RstFormat,
											  node.Patient,
											  node.TxInitiation,
											  node.TxRoom,
											  node.BAMS,
											  node.table,
											  node.gantry),
											  nodes.IES))
							
	node.PTTxPlan <- XML::xmlNode( "PTTxPlan", 
							  .children = list(node.Beam))

	node.BeamPlan <- XML::xmlNode( "Parameter", 
							  attrs = c(Name = "BeamPlan"),
							  .children = list(node.PTTxPlan))

	if(user.input$IES[[1]]$spill.intensity == 0){
    node.SfpParameters <- XML::xmlNode( "SfpParameters", 
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
	node.SfpParameters <- XML::xmlNode( "SfpParameters", 
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
	node.WorkflowContext <- XML::xmlNode( "WorkflowContext",
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
  node.PTTxPlanMd5 <- XML::xmlNode("PTTxPlanMd5",
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
					")\n     S. GREILICH, S. STEFANOWICZ, DKFZ -->\n",
					sep = "")
			 
	if(plan.format == 1){
	  XML::saveXML( node.Message, 
				 file = file.path(paste(path.save, "/", file.name.request, sep = "")),
				 prefix = prefix)

	  XML::saveXML( node.WorkflowContext, 
				 file = file.path(paste(path.save, "/", file.name.plan, sep = "")),
				 prefix = prefix)

		cat("\nSaved request file ",
			file.name.request,
			" and plan file ",
			file.name.plan,
			"\nto location ",
			path.save,
			"\n\nYou can directly copy this folder onto Rtt-Pt-Sim drive E:\n!!! IMPORTANT: If your plan does not work at HIT, esp. with only one ion species double check with accelerator team if MD% checksums for the control system have been updated after any kind of system maintainance. This is the most likely cause for HITXML plans to fail. !!!\nDone.\n",
			sep = "")
	}else{
    
	  # Save plan in the new format, according to Steffen's format (Yukihara, Oct. 2014)
	  XML::saveXML( node.PTTxPlanMd5, 
	           file = file.path(paste(path.save, "/", file.name.plan, sep = "")),
	           prefix = prefix)
    
#     # Save plan in simple new format, according to Julia's example (Yukihara, Oct. 2014)
# 	  saveXML( node.PTTxPlan, 
# 	           file = file.path(paste(path.save, "/", file.name.plan, sep = "")),
# 	           prefix = prefix)

cat("\nSaved plan file ",
			file.name.plan,
			"\nto location ",
			path.save,
			"\n\nYou can directly copy this folder onto Rtt-Pt-Sim drive E:\n!!! IMPORTANT: If your plan does not work at HIT, esp. with only one ion species double check with accelerator team if MD% checksums for the control system have been updated after any kind of system maintainance. This is the most likely cause for HITXML plans to fail. !!!\nDone.\n",
			sep = "")
	}
	
	if(!interactive.run.mode){
	  cat( "\n\nSave same plan with different name [y/n (default, press return)]?\n--> ",
			   sep = "")
	  answer <- readLines(con     = input.con,
	                      n       = 1)
	  if(answer != 'y'){
		  break
	  }	
	}else{
    break
	}
	# Get new experiment run name
	result          <- HX.prompt.user( variable.name = "name.exp.run",
									   prompt        = "Please give NEW experiment run name - unique, no whitespace",
									   df.defaults   = df.defaults)
	name.exp.run    <- result$value
	df.defaults     <- result$df.defaults	
}



# Start pdf viewer for results, Acrobat Reader will
# give problems here if alread open with same file
# Rather use a viewer as "evince" as default
# Check OS
if(!interactive.run.mode){
  current.OS <- "Linux"
  if (grepl("Windows", Sys.info()[[1]]) == TRUE){
      current.OS <- "Windows"
  }
  switch( match(current.OS, c("Linux", "Windows", "Apple")),
      system(paste("evince ./",path.save,"/",user.input$basic$name.exp.run, ".pdf", sep = "")),
  	system(paste("cmd /c start ", path.save,"\\",user.input$basic$name.exp.run, ".pdf", sep = "")),
      system("Tba"))
}
