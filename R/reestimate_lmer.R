#' Re-estimate model  
#' @param object ebp object 
#' Outputs lmerMod object 
#' 
reestimate_lmer <- function(object=object) {
  data <- object$model$data  
  formula <- object$model$call$fixed
  formula <- paste0(formula," + (1 | ",names(object$model$call$random[1]),")")
  if (!object$model$call$random[2]=="NULL") {
    formula <- paste0(formula," + (1 | ",names(object$model$call$random[2]),")")
  }
  args=list(formula=formula,data=data,weights=object$model$data[,object$framework$weight])
  
  lme <- do.call(lme4:::lmer,args)
  return(lme)
  }
