.onLoad <- function(lib, pkg){
 library.dynam("HITXML", pkg, lib)
 packageStartupMessage("HITXML package loaded.\n")
}

.onUnload <- function(libpath){
 try(library.dynam.unload("HITXML", libpath))
}

HITXML <- function(){
  source( system.file( 'exec', '1-Templates', 'HIT_XML.Create_Plan.R', package = 'HITXML' ) )
}