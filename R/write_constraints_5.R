write_constraints_5 <- function(variables = variables,
                                background.networks.list = background.networks.list,
                                tf_scores = tf_scores){
  
  constraints <- c()
  
  cell_types <- names(background.networks.list$background.networks)
  
  all_ligands <- background.networks.list$ligands.receptors$ligands
  
  # Constraints 5a
  cc1 <- c()
  cc2 <- c()
  for(ii in 1:length(all_ligands)){
    
    lr_var <- variables$var[which(variables$var_exp==paste0("LR:ligand ", all_ligands[ii]))]
    
    vv <- c()
    for(jj in 1:length(cell_types)){
      
      var <- variables$var[which(grepl(pattern = paste0(cell_types[jj], ":"), x = variables$var_exp, fixed = TRUE))]
      var_exp <- variables$var_exp[which(grepl(pattern = paste0(cell_types[jj], ":"), x = variables$var_exp, fixed = TRUE))]
      
      idx1 <- which(var_exp==paste0(cell_types[jj], ":node ", all_ligands[ii]))
      idx2 <- which(paste0(background.networks.list$background.networks[[jj]]$gene_source, 
                           "=",
                           background.networks.list$background.networks[[jj]]$gene_target) %in% 
                      paste0(tf_scores[[jj]]$tf, "=", all_ligands[ii]))
      if((length(idx1)==1) && (length(idx2)>0)){
        vv <- c(vv, var[idx1])
      }
      
    }
    
    cc1 <- c(cc1, paste0(length(cell_types), " ", lr_var, " - ", paste0(vv, collapse = " - "), " >= 0"))
    cc2 <- c(cc2, paste0(length(cell_types), " ", lr_var, " - ", paste0(vv, collapse = " - "), " <= ", length(cell_types)-1))
    
  }
  
  
  # Constraints 5b
  cc3 <- c()
  cc4 <- c()
  for(ii in 1:length(all_ligands)){
    
    for(jj in 1:length(cell_types)){
      
      var <- variables$var[which(grepl(pattern = paste0(cell_types[jj], ":"), x = variables$var_exp, fixed = TRUE))]
      var_exp <- variables$var_exp[which(grepl(pattern = paste0(cell_types[jj], ":"), x = variables$var_exp, fixed = TRUE))]
      
      intint <- paste0(cell_types[jj], ":interaction ", tf_scores[[jj]]$tf, "=", all_ligands[ii])
      
      idx1 <- which(var_exp==paste0(cell_types[jj], ":node ", all_ligands[ii]))
      idx2 <- which(var_exp %in% intint)
      
      
      if((length(idx1)==1) && (length(idx2) > 0)){
        cc3 <- c(cc3, paste0(var[idx1], " - ", paste0(var[idx2], collapse = " - "), " <= 0"))
        cc4 <- c(cc4, paste0(var[idx1], " - ", var[idx2], " >= 0"))
      }
      
    }
    
  }
  
  # Constraints 5c (a LR cannot happen if a ligand is not present in EC)
  cc5 <- c()
  for(ii in 1:length(all_ligands)){
    
    curr <- all_ligands[ii]
    lVar <- variables$var[which(variables$var_exp==paste0("LR:ligand ", curr))]
    ind <- which(grepl(pattern = paste0("LR:interaction ", curr, "="), 
                       x = variables$var_exp, fixed = TRUE))
    
    if(length(ind) > 0){
      
      cc5 <- c(cc5, paste0(variables$var[ind], " - ", lVar, " <= 0"))
      
    }
    
  }
  
  
  constraints <- c(cc1, cc2, cc3, cc4, cc5)
  
  return(constraints)

  
}