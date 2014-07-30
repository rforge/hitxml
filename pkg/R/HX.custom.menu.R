HX.custom.menu <- function(prompt, choices, default, input.con){
   cat("\n", prompt, "\n")
   for (i in 1:length(choices)){
      if(i == as.numeric(default)){
      	cat("   ", i, ": ", choices[i], " (default, press return to choose)\n", sep = "")
      }else{
      	cat("   ", i, ": ", choices[i], "\n", sep = "")
      }
   }
   cat("\n--> ")

  repeat{
    chosen <- readLines(con     = input.con,
                        n       = 1)
	if(chosen == ""){
		break
	}else{
		chosen <- as.numeric(chosen)
	}
    if(is.na(chosen)){
      chosen <- 0
    }
    if((chosen >= 1) & (chosen <= length(choices))){
      break
    }else{
      cat("Please enter a valid number or press return for default.\n: ")
    }
  }
  return(chosen)
}
