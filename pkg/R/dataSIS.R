#' @title dataSIS structure
#' 
#' @description Represents SIS (raster scanner information) data
#'
#' @author Steffen Greilich 
setClass( Class            = "dataSIS",
          slots            = c( projectile             = "character",
                                file.date              = "character",
                                beam.energies.MeV.u    = "numeric",
                                foci                   = "data.frame",
                                intensities            = "data.frame"),
          prototype        = list( projectile             = character(),
                                   beam.energies.MeV.u    = numeric(),
                                   foci                   = data.frame(beam.energy.MeV.u = numeric(),
                                                                       focus.no          = integer(),
                                                                       foces.FWHM.mm     = numeric()),
                                   intensities            = data.frame(beam.energy.MeV.u = numeric(),
                                                                       intensity.no      = integer(),
                                                                       intensity.N.s     = numeric())))

dataSIS <- function(file.name, sis.path){
  input		<-	scan( file.path(sis.path,  
                              file.name, 
                              fsep = .Platform$file.sep),
                 what = "character", 
                 sep = "\n")

  file.date     <-  read.item.character(input, "sistable")
  projectile	  <-	read.item.character(input, "projectile")

  sis           <- read.table( file.path(sis.path,  
                                         file.name, 
                                         fsep = .Platform$file.sep), 
                               skip = 5,
                               stringsAsFactors = FALSE)
  sis           <- sis[-nrow(sis),]    # delete last row
  names(sis)    <- as.character(1:ncol(sis))
  energy.row    <- grep("energy", sis[1,])+1
  focus.rows    <- seq( grep("focus", sis[1,])+1,
                        grep("intensity", sis[1,])-1,
                        by = 1)
  intensity.rows <- seq(grep("intensity", sis[1,])+1,
                        ncol(sis),
                        by = 1)
  
  
  foci        <- sis[,c(energy.row,focus.rows)]
  names(foci) <- c("beam.energy.MeV.u", as.character(1:length(focus.rows)))
  foci        <- tidyr::gather_(foci, "focus.no", "focus.FWHM.mm", names(foci)[-1])
  
  intensities        <- sis[,c(energy.row,intensity.rows)]
  names(intensities) <- c("beam.energy.MeV.u", as.character(1:length(intensity.rows)))
  intensities        <- tidyr::gather_(intensities, "intensity.no", "intensity.N.s", names(intensities)[-1])
  

  new("dataSIS",
      projectile          = projectile,
      file.date           = file.date,
      beam.energies.MeV.u = unlist(sis[,c(energy.row)]),
      foci                = foci,
      intensities         = intensities)
}

get.IES.index <- function(sis.data, beam.energy.MeV.u){
  which(abs(sis.data@beam.energies.MeV.u - beam.energy.MeV.u) == min(abs(sis.data@beam.energies.MeV.u - beam.energy.MeV.u)))
}

get.clostest.beam.energy.MeV.u <- function(sis.data, beam.energy.MeV.u){
  sis.data@beam.energies.MeV.u[get.IES.index(sis.data, beam.energy.MeV.u)]
}

get.foci <- function(sis.data, beam.energy.MeV.u){
  energy.MeV.u <- get.clostest.beam.energy.MeV.u(sis.data, beam.energy.MeV.u)
  ii           <- sis.data@foci$beam.energy.MeV.u == energy.MeV.u
  return(sis.data@foci[ii,])
}

get.focus.FWHM.mm <- function(sis.data, beam.energy.MeV.u, focus.no){
  ff <- get.foci(sis.data, beam.energy.MeV.u)
  ii <- ff$focus.no == focus.no
  return(ff$focus.FWHM.mm[ii])
}

get.intensities <- function(sis.data, beam.energy.MeV.u){
  energy.MeV.u <- get.clostest.beam.energy.MeV.u(sis.data, beam.energy.MeV.u)
  ii           <- sis.data@foci$beam.energy.MeV.u == energy.MeV.u
  return(sis.data@intensities[ii,])
}

get.intensity.N.s <- function(sis.data, beam.energy.MeV.u, intensity.no){
    ii <- get.intensities(sis.data, beam.energy.MeV.u)
    jj <- ii$intensity.no == intensity.no
    return(unlist(ii$intensity.N.s[jj,]))
}
