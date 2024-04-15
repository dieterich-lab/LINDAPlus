write_constraints_8 <- function(variables = variables, 
                                background.networks.list = background.networks.list){
  
  if(length(background.networks.list$ligands.receptors$perturbation.ligands) > 0){
    
    var_exp <- variables$var_exp[which(grepl(pattern = "PseudoCell:", x = variables$var_exp, fixed = TRUE))]
    vv <- variables$var[which(grepl(pattern = "PseudoCell:", x = variables$var_exp, fixed = TRUE))]
    
    ind <- which(grepl(pattern = "PseudoCell:dist ", x = var_exp, fixed = TRUE))
    var_exp <- var_exp[-ind]
    vv <- vv[-ind]
    
    constraints <- paste0(vv, " = 1")
    
  } else {
    
    constraints <- NULL
    
  }
  
  return(constraints)
  
}