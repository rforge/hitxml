.onLoad <- function(lib, pkg){
 library.dynam("HITXML", pkg, lib)
 packageStartupMessage("HITXML package loaded.\n")
}

.onUnload <- function(libpath){
 try(library.dynam.unload("HITXML", libpath))
}

HITXML <- function(){
  cat( system.file( 'exec', 'HIT_XML.R', package = 'HITXML' ) )
}