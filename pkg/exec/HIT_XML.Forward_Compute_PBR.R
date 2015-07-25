#!/usr/bin/Rscript

############################################
# Script to forward calculate HIT 
# physical beam records 
#
# S. Greilich, Apr. 2013
# DKFZ Heidelberg, s.greilich@dkfz.de
############################################

# Clear workspace, set version
rm(list =ls())

# load package
library(HITXML)
program.version <- packageDescription("HITXML")$Version
program.date    <- packageDescription("HITXML")$Date

# Say hello to user
welcome.message <- paste( "\n",
                          "#####################################################################\n",
                          "#####################################################################\n",
                          "## THIS IS HIT_XML_PBR\n## S. Greilich, DKFZ Heidelberg\n## Version ", 
                          program.version, 
                          " as of ", 
                          program.date,
                          "\n",
                          "#####################################################################\n",
                          sep = "")

cat(welcome.message)

# Script run interactively or non-interactively, set connection accordingly
interactive.run.mode <- interactive()
if(!interactive.run.mode){
  input.con <-file("stdin")
}else{
  input.con <-stdin()
}

resolution.mm          <- as.numeric(HX.prompt.user( prompt        = "Resolution for optimizing homogeneity - in mm, depending e.g. on your detector resolution (if unsure \'1\' is a good guess",
                                                                      default   = 1,
                                                                      input.con = input.con))

# Select beam record
files              <- list.files( path = ".",
                                  pattern = "xml$")

repeat{
  
  repeat{
    PBR.file           <- as.numeric(HX.prompt.user( prompt        = "Please choose PBR file",
                                       choices       = c(files, "None, quit"),
                                       default   = 1,
                                       input.con = input.con)) 
    if(PBR.file == length(files) + 1){
      stop("Quitting...")
    }
    file.name          <- files[PBR.file]
    
    cat("Chosen:", file.name)
    
    
    if(!interactive.run.mode){
      file.out <- paste(file.name, ".pdf", sep = "")
      pdf(file = file.out)
    }
  
    # DEBUG:
    # file.name <- "D:/workspaces/R/Forward/PhysicalBeamRecord_TCU1_20130222020244.xml"
    df <- HX.read.PBR(file.name)

    if(min(df$focus.X.FWHM.mm) <= 0 || min(df$focus.Y.FWHM.mm) <= 0){
      cat("\n!! PBR contains no focus information (TCU3?), skipping...\n\n")
      dev.off()
    }else{
      break
    }
  }
  
  n.IES <- unique(df$IES)
  
  # Find max field
  field.mm    <- HX.getMinAndMax(beam.spot.grid = df)
  
  
  # Add all IESs if applicable
  if(n.IES > 1){
    stop("Multi-IES not yet implemented")
    # Should be tapply
    #for(i in 1:n.IES){
      # i <- 1
     # tmp <- HX.compute.field( beam.spot.grid = df[df$IES == i,],
      #                         field.mm       = field.mm,
      #                         resolution.mm  = resolution.mm)
    #}
  }
  beam.spot.grid <- df[df$IES == 1,]
  beam.spot.grid$spot.no <- 1:nrow(beam.spot.grid)
  
  m   <- HX.compute.field( beam.spot.grid = beam.spot.grid,
                           field.mm       = field.mm,
                           resolution.mm  = resolution.mm)
  
  title.doc <- paste(welcome.message,
                     "\n\nHIT_XML_PBR forward calculation\n\n",
                     "PBR name: ",
                     file.name,
                     "\n",
                     sep = "")
  
  clan.txtplot(title.doc)
  
  # Plot field
  require(lattice)
  
  # Plot spot positions
  custom.panel    <- function(spot.no = spot.no, x = x, y = y, ...){
    panel.xyplot(x = x, y = y, ...)
    for(i in 1:length(spot.no)){
      panel.text( spot.no[i],
                  x = x[i],
                  y = y[i],
                  cex = 10 / sqrt(nrow(beam.spot.grid)))
    }
  }
  
  plot(
    xyplot( y.mm ~ x.mm,
            beam.spot.grid,
            type = 'p',
            pch  = 3,
            cex  = 10 / sqrt(nrow(beam.spot.grid)),
            spot.no = beam.spot.grid$spot.no,
            panel = custom.panel,
            xlim = c(min(m[,1]), max(m[,1])),
            ylim = c(min(m[,2]), max(m[,2])),
            main = paste('entire field (IES no. ', n.IES, ')', sep = ""), 
            sub  = 'beam spot positions and numbers',
            aspect = 1)
  )		
  # Plot full field
  
  y.line.mm       <- min(abs(m[,2]))
  custom.panel    <- function(...){
    panel.levelplot(...)
    panel.xyplot( beam.spot.grid$x.mm, 
                  beam.spot.grid$y.mm,
                  type = 'p',
                  pch  = 3)
    panel.abline(h = y.line.mm, lty = 2)
  }
  
  plot(
    levelplot(m[,3]*100 ~ m[,1]*m[,2],
              panel     = custom.panel,
              cex       = 20 / sqrt(nrow(beam.spot.grid)),
              xlab      = 'x / mm',
              ylab      = 'y / mm',
              main      = paste('entire field (IES no. ', n.IES, ')', sep = ""),
              sub       = 'fluence / (1/cm2)',
              aspect    = 1)
  )
  
  # Plot cross section
  plot(
    xyplot(m[,3]*100 ~ m[,1],
           type   = 's',
           subset = abs(m[,2]) == y.line.mm,
           xlab   = 'x / mm',
           ylab   = 'fluence / (1/cm2)',
           main   = paste('cross section (IES no. ', n.IES, ')', sep = ""))
  )
  
  # Plot homogenous part, with relative deviation from requested fluence
  m               <- HX.compute.field( beam.spot.grid = beam.spot.grid, 
                                       field.mm       = c( min(beam.spot.grid$x.mm) + max(beam.spot.grid$focus.X.FWHM.mm) * 2, 
                                                           min(beam.spot.grid$y.mm) + max(beam.spot.grid$focus.Y.FWHM.mm) * 2, 
                                                           max(beam.spot.grid$x.mm) - max(beam.spot.grid$focus.X.FWHM.mm) * 2, 
                                                           max(beam.spot.grid$y.mm) - max(beam.spot.grid$focus.Y.FWHM.mm) * 2),
                                       resolution.mm  = resolution.mm)
  
  fluence.cm2 <- mean(m[,3])
  plot(
    levelplot( (m[,3] - fluence.cm2) /fluence.cm2 * 100 ~ m[,1]*m[,2],
               panel = custom.panel,
               xlab  = 'x / mm',
               ylab  = 'y / mm',
               main      = paste('homogenous central area (IES no. ', n.IES, ')', sep = ""),
               sub       = 'deviation from mean fluence / %'),
  )
  
  
  if(!interactive.run.mode){
    cat( "\n\nProcess next PBR [y/n (default, press return)]?\n--> ",
         sep = "")
    answer <- readLines(con     = input.con,
                        n       = 1)
    if(answer != 'y'){
      break
    }	
  }else{
    break
  }
  dev.off()
}# Repeat



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
          system(paste("evince ", file.out, ".pdf", sep = "")),
          system(paste("cmd /c start ", file.out, ".pdf", sep = "")),
          system("Tba"))
}
