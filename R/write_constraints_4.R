write_constraints_4 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  constraints <- c()
  
  cell_types <- names(background.networks.list$background.networks)
  
  # IC Part
  for(ii in 1:length(cell_types)){
    
    var <- variables$var[which(grepl(pattern = paste0(cell_types[ii], ":"), x = variables$var_exp, fixed = TRUE))]
    var_exp <- variables$var_exp[which(grepl(pattern = paste0(cell_types[ii], ":"), x = variables$var_exp, fixed = TRUE))]
    var_exp <- gsub(pattern = paste0(cell_types[ii], ":"), replacement = "", x = var_exp, fixed = TRUE)
    
    cc1 <- c()
    cc2 <- c()
    idx <- which(grepl(pattern = "interaction ", x = var_exp, fixed = TRUE))
    for(jj in 1:length(idx)){
      
      curr_int <- gsub(pattern = "interaction ", replacement = "", x = var_exp[idx[jj]], fixed = TRUE)
      
      ind <- which(sapply(strsplit(var_exp, split = " ", fixed = TRUE), "[", 4) == curr_int)
      
      if(length(ind) > 0){
        
        cc1 <- c(cc1, paste0(var[idx[jj]], " - ", var[ind], " >= 0"))
        cc2 <- c(cc2, paste0(var[idx[jj]], " - ", paste0(var[ind], collapse = " - "), " <= 0"))
        
      }
      
    }
    
    constraints <- c(constraints, cc1, cc2)
    
  }
  
  #EC Part
  var <- variables$var[which(grepl(pattern = "LR:", x = variables$var_exp, fixed = TRUE))]
  var_exp <- variables$var_exp[which(grepl(pattern = "LR:", x = variables$var_exp, fixed = TRUE))]
  var_exp <- gsub(pattern = "LR:", replacement = "", x = var_exp, fixed = TRUE)
  reac <- var_exp[which(grepl(pattern = "reaction ", x = var_exp, fixed = TRUE))]
  reac <- gsub(pattern = "reaction ", replacement = "", x = reac, fixed = TRUE)
  ppi <- paste0(sapply(strsplit(x = sapply(strsplit(x = reac, split = "=", fixed = TRUE), "[", 1), split = "_", fixed = TRUE), "[", 2),
                "=",
                sapply(strsplit(x = sapply(strsplit(x = reac, split = "=", fixed = TRUE), "[", 2), split = "_", fixed = TRUE), "[", 2))
  
  cc1 <- c()
  cc2 <- c()
  idx <- which(grepl(pattern = "interaction ", x = var_exp, fixed = TRUE))
  for(jj in 1:length(idx)){
    
    curr_int <- gsub(pattern = "interaction ", replacement = "", x = var_exp[idx[jj]], fixed = TRUE)
    
    ind <- which(ppi==curr_int)
    
    if(length(ind) > 0){
      
      cc1 <- c(cc1, paste0(var[idx[jj]], " - ", var[which(var_exp%in%paste0("reaction ", reac[ind]))], " >= 0"))
      cc2 <- c(cc2, paste0(var[idx[jj]], " - ", paste0(var[which(var_exp%in%paste0("reaction ", reac[ind]))], collapse = " - "), " <= 0"))
      
    }
    
  }
  
  constraints <- c(constraints, cc1, cc2)
  
  return(constraints)
  
}