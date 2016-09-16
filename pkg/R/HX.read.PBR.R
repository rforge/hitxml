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
  
  df <- NULL
  for (i in 1:n.IES) {
    cur.IES <- df.IES$IES[i]
    x       <- as.numeric(xpathApply(doc, 
                                     paste0("//IES[@number=\"", 
                                            cur.IES, 
                                            "\"]/Voxel"), 
                                     xmlGetAttr,
                                     "x"))
    y       <- as.numeric(xpathApply(doc, 
                                     paste0("//IES[@number=\"", 
                                            cur.IES, 
                                            "\"]/Voxel"), 
                                     xmlGetAttr,
                                     "y"))
    N       <- as.numeric(xpathApply(doc, 
                                     paste0("//IES[@number=\"", 
                                            cur.IES, 
                                            "\"]/Voxel"), 
                                     xmlGetAttr,
                                     "particles"))
    df      <- rbind(df, data.frame(IES         = cur.IES,
                                    x.mm        = x,
                                    y.mm        = y,
                                    N.particles = N))
  }
  df <- merge(df, df.IES, by = "IES")
  return(df)
}
