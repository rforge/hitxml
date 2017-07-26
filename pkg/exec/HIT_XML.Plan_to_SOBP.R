rm(list = ls())

library(XML)

doc <- XML::xmlTreeParse("bio_c_1Gy_8mm_35mm_1x25_schmal.xml", useInternalNodes = TRUE)

IESs <- XML::getNodeSet(doc, path = "//IES")
ll   <- lapply(IESs,
       function(x){
         # x <- IESs[[1]]
         energy.GeV <- as.numeric(XML::xmlGetAttr(x, name = "energy"))/1000.0
         focus.cm   <- as.numeric(XML::xmlGetAttr(x, name = "focus"))/10.0
         pos.x.cm   <- sapply(XML::getNodeSet(x, path = "//Voxel"),
                              function(y){
                                as.numeric(XML::xmlGetAttr(y, name = "x")) / 10.0
                              })
         pos.y.cm   <- sapply(XML::getNodeSet(x, path = "//Voxel"),
                              function(y){
                                as.numeric(XML::xmlGetAttr(y, name = "y")) / 10.0
                              })
         particles  <- sapply(XML::getNodeSet(x, path = "//Voxel"),
                              function(y){
                                as.numeric(XML::xmlGetAttr(y, name = "particles")) / 10.0
                              })
         return(data.frame(energy.GeV = energy.GeV,
                           pos.x.cm = pos.x.cm,
                           pos.y.cm = pos.y.cm,
                           focus.cm = focus.cm,
                           particles = particles))
       })

df <- do.call("rbind", ll)

write.table(df,
            file = "SOBP_from_Plan.dat",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol       = "\r\n" )
