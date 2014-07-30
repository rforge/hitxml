HX.prompt.user <- function(prompt, choices = NULL, default, input.con){

    # Simple prompt or menu?
    if(is.null(choices)){
	    cat( "\n",
		     prompt,
			 " (return for \"", 
		     default, 
			 "\")\n--> ",
			 sep = "")
		answer <- readLines(con     = input.con,
		                    n       = 1)
    }else{
    	answer <- HX.custom.menu(prompt        = prompt, 
					             choices       = choices, 
					             default       = default,
                       input.con     = input.con)
    }

  if(answer == ""){
	    value                  <- default
	}else{
		  value                  <- answer
	}
    
  return(value)
}
