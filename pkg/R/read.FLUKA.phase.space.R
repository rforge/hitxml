#' @title Import of phase space files scored with FLUKA
#'
#' @description Reads in phase space data scored with the mgdraw routine supplied with
#' this package, can be from multiple runs and multiple boundaries
#'
#' @param path path containing the phase space data
#' @param pattern file pattern
#' @param particles list of FLUKA particle names that will be kept, if NULL, all particles will be kept
#' @param remove.backscattered if true, particles travelling upstream are removed
#' @param Z.max optional upper limit for particle charge
#'
#' @return data frame
#'
#' @author Steffen Greilich
read.FLUKA.phase.space <- function(path,
                                   pattern = "*_phase_space",
                                   particles = c("PROTON", "DEUTERON", "TRITON",
                                                 "3-HELIUM", "4-HELIUM",
                                                 "HEAVYION"),
                                   remove.backscattered = TRUE,
                                   region.name.prefix = "TARGET",
                                   Z.max) {
  file.list <- list.files(path,
                          pattern = pattern,
                          full.names = TRUE)
  
  data.list <- lapply(file.list,
                      function(x) {
                        cat("Reading", x, "...\n")
                        read.table(x,
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
                      })
  
  df.raw        <- do.call("rbind", data.list)
  df.raw$travel <-
    paste0("From: ", df.raw$OLDREG, " / to: ", df.raw$NEWREG)
  
  # There are actually some particles being backscattered...
  # table(df.raw$travel)
  # Remove them (in a pretty stupid way)
  if (remove.backscattered) {
    ii <-
      as.numeric(gsub(region.name.prefix, "", df.raw$OLDREG)) <  as.numeric(gsub(region.name.prefix, "", df.raw$NEWREG))
    cat(
      "Removing entries for backscattered particles (", sum(!ii), "of", nrow(df.raw), ")\n"
    )
    df.raw <- df.raw[ii,]
  }
  
  if (!is.null(particles)) {
    # Filter out unwanted particles
    ii     <- df.raw$PRNAME %in% particles
    gg     <- table(df.raw$PRNAME[!ii])
    cat("Removing particles other than", paste(particles, collapse = ", "), ",\n")
    cat("   which are", sum(!ii), "of", nrow(df.raw), "entries\n")
    cat("Table of particles removed:\n")
    print(gg)
    df.raw <- df.raw[ii,]
  }
  
  if (!missing(Z.max)) {
    ii     <- df.raw$ICHRGE <= Z.max
    cat("Removing particles with charge higher than", Z.max, "(", sum(!ii), "of", nrow(df.raw), "entries)\n")
    df.raw <- df.raw[ii,]
  }
  
  return(df.raw)
}
