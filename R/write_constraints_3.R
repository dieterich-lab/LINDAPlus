write_constraints_3 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  
  constraints <- c()
  
  cell_types <- names(background.networks.list$background.networks)
  
  for(ii in 1:length(cell_types)){
    
    var <- variables$var[which(grepl(pattern = paste0(cell_types[ii], ":"), x = variables$var_exp, fixed = TRUE))]
    var_exp <- variables$var_exp[which(grepl(pattern = paste0(cell_types[ii], ":"), x = variables$var_exp, fixed = TRUE))]
    var_exp <- gsub(pattern = paste0(cell_types[ii], ":"), replacement = "", x = var_exp, fixed = TRUE)
    
    # Write the out-going constraints
    cc1 <- c()
    idx <- which(grepl(pattern = "node ", x = var_exp, fixed = TRUE))
    for(jj in 1:length(idx)){
      
      ind <- which(grepl(pattern = paste0("interaction ", gsub(pattern = "node ", replacement = "", fixed = TRUE, x = var_exp[idx[jj]]), "="), 
                         x = var_exp, fixed = TRUE))
      
      if(length(ind) > 0){
        
        cc1 <- c(cc1, paste0(paste0(var[ind], collapse = " + "), " - ", var[idx[jj]], " >= 0"))
        
      }
      
    }
    
    # Write the incoming constraints
    cc2 <- c()
    idx <- which(var_exp%in%paste0("node ", setdiff(unique(c(background.networks.list$background.networks[[ii]]$gene_source,
                                                             background.networks.list$background.networks[[ii]]$gene_target)),
                                                    background.networks.list$ligands.receptors$receptors)))
    for(jj in 1:length(idx)){
      
      ind1 <- which(grepl(pattern = "interaction ", x = var_exp, fixed = TRUE))
      ind2 <- which(sapply(strsplit(x = var_exp, split = "=", fixed = TRUE), "[", 2)==gsub(pattern = "node ", replacement = "", x = var_exp[idx[jj]], fixed = TRUE))
      
      ind <- intersect(x = ind1, y = ind2)
      
      if(length(ind) > 0){
        
        cc2 <- c(cc2, paste0(paste0(var[ind], collapse = " + "), " - ", var[idx[jj]], " >= 0"))
        
      }
      
    }
    
    constraints <- c(constraints, cc1, cc2)
    
  }
  
  
  return(constraints)
  
}