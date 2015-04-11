HX.read.PBR <- function (file.name, IES.only = FALSE) 
{
  doc <- xmlTreeParse(file.name, useInternalNodes = TRUE)
  n.IES <- length(getNodeSet(doc = doc, path = "//IES"))
  df.IES <- data.frame(IES = numeric(n.IES), energy.MeV.u = numeric(n.IES), 
                       focus.X.FWHM.mm = numeric(n.IES), focus.Y.FWHM.mm = numeric(n.IES))
  df.IES$IES <- as.numeric(xpathApply(doc, "//IES", xmlGetAttr, 
                                      "number"))
  df.IES$energy.MeV.u <- as.numeric(xpathApply(doc, "//IES", 
                                               xmlGetAttr, "energy"))
  if (length(getNodeSet(doc = doc, path = "//PTTxRecord")) == 1) {
    df.IES$focus.X.FWHM.mm <- as.numeric(xpathApply(doc, 
                                                    "//IES", xmlGetAttr, "focusX"))
    df.IES$focus.Y.FWHM.mm <- as.numeric(xpathApply(doc, 
                                                    "//IES", xmlGetAttr, "focusY"))
  } else {
    df.IES$focus.X.FWHM.mm <- as.numeric(xpathApply(doc, 
                                                    "//IES", xmlGetAttr, "focus"))
    df.IES$focus.Y.FWHM.mm <- df.IES$focus.X.FWHM.mm
  }
  if (IES.only == TRUE) {
    return(df.IES)
  }
  n.voxel <- length(getNodeSet(doc = doc, path = "//Voxel"))
  df <- data.frame(IES = numeric(n.voxel), x.mm = numeric(n.voxel), 
                   y.mm = numeric(n.voxel), N.particles = numeric(n.voxel))
  cur.idx <- 1
  for (i in 1:n.IES) {
    cur.IES <- df.IES$IES[i]
    tmp.x <- as.numeric(xpathApply(doc, paste("//IES[@number=\"", 
                                              cur.IES, "\"]/Voxel", sep = ""), xmlGetAttr, "x"))
    idx <- cur.idx:(cur.idx + length(tmp.x) - 1)
    df$IES[idx] <- rep(cur.IES, length(idx))
    df$x.mm[idx] <- tmp.x
    df$y.mm[idx] <- as.numeric(xpathApply(doc, paste("//IES[@number=\"", 
                                                     cur.IES, "\"]/Voxel", sep = ""), xmlGetAttr, "y"))
    df$N.particles[idx] <- as.numeric(xpathApply(doc, paste("//IES[@number=\"", 
                                                            cur.IES, "\"]/Voxel", sep = ""), xmlGetAttr, "particles"))
  }
  df <- merge(df, df.IES, by = "IES")
  return(df)
}
