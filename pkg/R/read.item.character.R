read.item.character  <-  function(data, code){
  # returns string w/o leading or trailing whitespace
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  line  <-	data[grep(code, data)]
  if (!is.null(line)){
    return(trim(substring( line, regexpr(" ", line) + 1, nchar(line))))
  }else{
    return(0.0)
  }
}
