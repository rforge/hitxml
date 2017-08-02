rm(list = ls())

library(XML)
library(ggplot2)
library(data.table)

doc <- XML::xmlTreeParse("bio_c_1Gy_8mm_35mm_1x25_schmal.xml", useInternalNodes = TRUE)

IESs <- XML::getNodeSet(doc, path = "//IES")
ll   <- lapply(IESs,
       function(x){
         # x <- IESs[[1]]
         energy.GeV <- as.numeric(XML::xmlGetAttr(x, name = "energy"))/1000.0
         focus.cm   <- as.numeric(XML::xmlGetAttr(x, name = "focus"))/10.0
         pos.x.cm   <- sapply(XML::xmlElementsByTagName(x, name = "Voxel"),
                              function(y){
                                as.numeric(XML::xmlGetAttr(y, name = "x")) / 10.0
                              })
         pos.y.cm   <- sapply(XML::xmlElementsByTagName(x, name = "Voxel"),
                              function(y){
                                as.numeric(XML::xmlGetAttr(y, name = "y")) / 10.0
                              })
         particles  <- sapply(XML::xmlElementsByTagName(x, name = "Voxel"),
                              function(y){
                                as.numeric(XML::xmlGetAttr(y, name = "particles"))
                              })
         return(data.frame(energy.GeV = energy.GeV,
                           pos.x.cm = pos.x.cm,
                           pos.y.cm = pos.y.cm,
                           focus.cm = focus.cm,
                           particles = particles))
       })

df <- as.data.table(do.call("rbind", ll))

tapply(df$particles, df$energy.GeV, unique)

df2 <- df[,
       .(energy.GeV = unique(energy.GeV),
         pos.x.cm  = 0.0,
         pos.y.cm  = 0.0,
         focus.cm  = unique(focus.cm),
         total.particles = sum(particles)),
       by = .(energy.GeV)]

ggplot(df, aes(x = pos.x.cm, y = pos.y.cm)) +
  theme_bw() +
  geom_point() +
  facet_wrap("energy.GeV")


write.table(df,
            file = "SOBP_from_Plan.dat",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol       = "\r\n" )

write.table(df2,
            file = "SOBP_from_Plan_simple.dat",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol       = "\r\n" )