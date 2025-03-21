write_constraints_5a <- function(variables = variables,
                                 background.networks.list = background.networks.list,
                                 tf.input = tf.input,
                                 constraits.parallel.writing = FALSE) {
  
  constraints <- c()
  cell_types <- names(background.networks.list$background.networks)
  all_ligands <- setdiff(background.networks.list$ligands.receptors$ligands, "PSEUDOLIGAND")
  
  # Precompute expressions for faster access
  lr_ligand_exprs <- paste0("LR:ligand ", all_ligands)
  lr_ligand_vars <- variables$var[match(lr_ligand_exprs, variables$var_exp)]
  
  process_ligand_constraints <- function(ii) {
    ligand <- all_ligands[ii]
    vv <- numeric(length(cell_types))
    
    for (jj in seq_along(cell_types)) {
      cell <- cell_types[jj]
      cell_prefix <- paste0(cell, ":")
      
      cell_var_idx <- grepl(cell_prefix, variables$var_exp)
      cell_vars <- variables$var[cell_var_idx]
      cell_var_exp <- sub(cell_prefix, "", variables$var_exp[cell_var_idx])
      
      idx1 <- match(paste0("node ", ligand), cell_var_exp)
      interactions <- paste0(tf.input[[jj]]$tf, "=", ligand)
      idx2 <- which(cell_var_exp %in% paste0("interaction ", interactions))
      
      if (!is.na(idx1) && length(idx2) > 0) {
        vv[jj] <- cell_vars[idx1]
      }
    }
    
    vv <- vv[vv != 0]  # Remove zero entries
    if (length(vv) > 0) {
      cc1 <- paste0(length(vv), " ", lr_ligand_vars[ii], " - ", paste(vv, collapse = " - "), " >= 0")
      cc2 <- paste0(length(vv), " ", lr_ligand_vars[ii], " - ", paste(vv, collapse = " - "), " <= ", length(vv) - 1)
      return(list(cc1 = cc1, cc2 = cc2))
    }
    return(NULL)
  }
  
  if (constraits.parallel.writing) {
    ligand_constraints <- mclapply(seq_along(all_ligands), process_ligand_constraints, mc.cores = length(cell_types))
  } else {
    ligand_constraints <- lapply(seq_along(all_ligands), process_ligand_constraints)
  }
  
  cc1 <- unlist(lapply(ligand_constraints, "[[", "cc1"))
  cc2 <- unlist(lapply(ligand_constraints, "[[", "cc2"))
  
  # Constraints 5c
  process_interaction_constraints <- function(ii) {
    ligand <- all_ligands[ii]
    ligand_var <- variables$var[match(paste0("LR:ligand ", ligand), variables$var_exp)]
    interaction_exprs <- paste0("LR:interaction ", ligand, "=")
    interaction_idxs <- which(grepl(interaction_exprs, variables$var_exp))
    
    if (length(interaction_idxs) > 0) {
      return(c(paste0(variables$var[interaction_idxs], " - ", ligand_var, " <= 0"),
               paste0(ligand_var, " - ", paste0(variables$var[interaction_idxs], collapse = " - "), " <= 0")))
    }
    return(NULL)
  }
  
  if (constraits.parallel.writing) {
    interaction_constraints <- mclapply(seq_along(all_ligands), process_interaction_constraints, mc.cores = length(cell_types))
  } else {
    interaction_constraints <- lapply(seq_along(all_ligands), process_interaction_constraints)
  }
  
  cc5 <- unlist(Filter(Negate(is.null), interaction_constraints))
  
  # Combine all constraints
  constraints <- c(cc1, cc2, cc5)
  
  return(constraints)
}
