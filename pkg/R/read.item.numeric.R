read.item.numeric  <-  function(data, code){
  line  <-	data[grep(code, data)]
  if (!is.null(line)){
    return(as.numeric(substring( line, regexpr(" ", line) + 1, nchar(line))))
  }else{
    return(0.0)
  }
}
