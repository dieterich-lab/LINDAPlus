write_constraints_5 <- function(variables = variables,
                                background.networks.list = background.networks.list,
                                tf.scores = tf.scores) {
  
  constraints <- c()
  
  cell_types <- names(background.networks.list$background.networks)
  all_ligands <- setdiff(background.networks.list$ligands.receptors$ligands, "PSEUDOLIGAND")
  
  # Precompute expressions for faster access
  lr_ligand_exprs <- paste0("LR:ligand ", all_ligands)
  lr_ligand_vars <- variables$var[match(lr_ligand_exprs, variables$var_exp)]
  
  # Constraints 5a
  cc1 <- cc2 <- vector("list", length(all_ligands))
  
  for (ii in seq_along(all_ligands)) {
    ligand <- all_ligands[ii]
    vv <- numeric(length(cell_types))
    
    for (jj in seq_along(cell_types)) {
      cell <- cell_types[jj]
      cell_prefix <- paste0(cell, ":")
      
      cell_var_idx <- grepl(cell_prefix, variables$var_exp)
      cell_vars <- variables$var[cell_var_idx]
      cell_var_exp <- sub(cell_prefix, "", variables$var_exp[cell_var_idx])
      
      idx1 <- match(paste0("node ", ligand), cell_var_exp)
      interactions <- paste0(tf.scores[[jj]]$tf, "=", ligand)
      idx2 <- which(cell_var_exp %in% paste0("interaction ", interactions))
      
      if (!is.na(idx1) && length(idx2) > 0) {
        vv[jj] <- cell_vars[idx1]
      }
    }
    
    vv <- vv[vv != 0]  # Remove zero entries
    if (length(vv) > 0) {
      cc1[[ii]] <- paste0(length(cell_types), " ", lr_ligand_vars[ii], " - ", paste(vv, collapse = " - "), " >= 0")
      cc2[[ii]] <- paste0(length(cell_types), " ", lr_ligand_vars[ii], " - ", paste(vv, collapse = " - "), " <= ", length(cell_types) - 1)
    }
  }
  
  # Flatten lists and filter NULLs
  cc1 <- unlist(Filter(Negate(is.null), cc1))
  cc2 <- unlist(Filter(Negate(is.null), cc2))
  
  # Constraints 5b and 5c - these could be similarly vectorized and optimized
  # Placeholder for actual optimized logic...
  
  # Constraints 5b
  cc3 <- c()
  cc4 <- c()
  for(ii in seq_along(all_ligands)){
    
    if(all_ligands[ii] != "PSEUDOLIGAND"){
      
      for(jj in seq_along(cell_types)){
        
        var <- variables$var[which(grepl(pattern = paste0(cell_types[jj], ":"), x = variables$var_exp, fixed = TRUE))]
        var_exp <- variables$var_exp[which(grepl(pattern = paste0(cell_types[jj], ":"), x = variables$var_exp, fixed = TRUE))]
        
        intint <- paste0(cell_types[jj], ":interaction ", tf.scores[[jj]]$tf, "=", all_ligands[ii])
        
        idx1 <- which(var_exp==paste0(cell_types[jj], ":node ", all_ligands[ii]))
        idx2 <- which(var_exp %in% intint)
        
        
        if((length(idx1)==1) && (length(idx2) > 0)){
          cc3 <- c(cc3, paste0(var[idx1], " - ", paste0(var[idx2], collapse = " - "), " <= 0"))
          cc4 <- c(cc4, paste0(var[idx1], " - ", var[idx2], " >= 0"))
        }
        
      }
      
    }
    
  }
  
  # Constraints 5c
  cc5 <- vector("list", length(all_ligands))
  for (ii in seq_along(all_ligands)) {
    ligand <- all_ligands[ii]
    ligand_var <- variables$var[match(paste0("LR:ligand ", ligand), variables$var_exp)]
    interaction_exprs <- paste0("LR:interaction ", ligand, "=")
    interaction_idxs <- which(grepl(interaction_exprs, variables$var_exp))
    
    if (length(interaction_idxs) > 0) {
      cc5[[ii]] <- paste0(variables$var[interaction_idxs], " - ", ligand_var, " <= 0")
    }
  }
  
  cc5 <- unlist(Filter(Negate(is.null), cc5))
  
  # Combine all constraints
  constraints <- c(cc1, cc2, cc3, cc4, cc5)
  
  return(constraints)
  
  constraints <- c(cc1, cc2)  # Include other constraints when optimized
  
  return(constraints)
}
