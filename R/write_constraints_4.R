write_constraints_4 <- function(variables = variables,
                                background.networks.list = background.networks.list,
                                constraits.parallel.writing = FALSE) {
  
  constraints <- vector(mode = "character")
  cell_types <- names(background.networks.list$background.networks)
  
  # Vectorized operations for the 'LR:' pattern matching and replacement
  lr_mask <- grepl("LR:", variables$var_exp, fixed = TRUE)
  var_lr <- variables$var[lr_mask]
  var_exp_lr <- sub("LR:", "", variables$var_exp[lr_mask], fixed = TRUE)
  reac <- var_exp_lr[grepl("reaction ", var_exp_lr, fixed = TRUE)]
  reac <- sub("reaction ", "", reac, fixed = TRUE)
  reac_splits <- strsplit(reac, "=", fixed = TRUE)
  ppi <- sapply(reac_splits, function(x) paste(sub("_.*", "", x[1]), sub(".*_", "", x[2]), sep = "="))
  
  process_cell_type <- function(cell) {
    local_constraints <- character()
    
    cell_mask <- grepl(paste0(cell, ":"), variables$var_exp, fixed = TRUE)
    var <- variables$var[cell_mask]
    var_exp <- sub(paste0(cell, ":"), "", variables$var_exp[cell_mask], fixed = TRUE)
    interaction_mask <- grepl("interaction ", var_exp, fixed = TRUE)
    interactions <- sub("interaction ", "", var_exp[interaction_mask], fixed = TRUE)
    interaction_vars <- var[interaction_mask]
    idx <- which(interaction_mask)
    
    if (length(idx) > 0) {
      interaction_details <- strsplit(var_exp, " ", fixed = TRUE)
      curr_ints <- sapply(interaction_details, function(x) x[4])
      cc1 <- cc2 <- character()
      
      for (jj in seq_along(idx)) {
        ind <- which(curr_ints == interactions[jj])
        if (length(ind) > 0) {
          cc1 <- c(cc1, paste0(var[idx[jj]], " - ", var[ind], " >= 0"))
          cc2 <- c(cc2, paste0(var[idx[jj]], " - ", paste0(var[ind], collapse = " - "), " <= 0"))
        }
      }
      
      local_constraints <- c(local_constraints, cc1, cc2)
    }
    
    var2rem <- variables$var[grepl(paste0(cell, ":reaction PSEUDODOMAINLR="), variables$var_exp, fixed = TRUE)]
    if (length(var2rem) > 0) {
      remrem <- unique(unlist(lapply(var2rem, function(x) which(grepl(x, local_constraints, fixed = TRUE)))))
      if (length(remrem) > 0) {
        local_constraints <- local_constraints[-remrem]
      }
    }
    
    return(local_constraints)
  }
  
  if (constraits.parallel.writing) {
    constraints_list <- mclapply(cell_types, process_cell_type, mc.cores = length(cell_types))
  } else {
    constraints_list <- lapply(cell_types, process_cell_type)
  }
  
  constraints <- unique(unlist(constraints_list))
  
  # EC Part
  var <- variables$var[grepl("LR:", variables$var_exp, fixed = TRUE)]
  var_exp <- sub("LR:", "", variables$var_exp[grepl("LR:", variables$var_exp, fixed = TRUE)], fixed = TRUE)
  reac <- sub("reaction ", "", var_exp[grepl("reaction ", var_exp, fixed = TRUE)], fixed = TRUE)
  ppi <- paste0(
    sapply(strsplit(sapply(strsplit(reac, "=", fixed = TRUE), "[", 1), "_", fixed = TRUE), "[", 2),
    "=",
    sapply(strsplit(sapply(strsplit(reac, "=", fixed = TRUE), "[", 2), "_", fixed = TRUE), "[", 2)
  )
  
  cc1 <- cc2 <- character()
  idx <- which(grepl("interaction ", var_exp, fixed = TRUE))
  for (jj in seq_along(idx)) {
    curr_int <- sub("interaction ", "", var_exp[idx[jj]], fixed = TRUE)
    ind <- which(ppi == curr_int)
    
    if (length(ind) > 0) {
      cc1 <- c(cc1, paste0(var[idx[jj]], " - ", var[which(var_exp %in% paste0("reaction ", reac[ind]))], " >= 0"))
      cc2 <- c(cc2, paste0(var[idx[jj]], " - ", paste0(var[which(var_exp %in% paste0("reaction ", reac[ind]))], collapse = " - "), " <= 0"))
    }
  }
  
  constraints <- c(constraints, cc1, cc2)
  
  return(constraints)
}
