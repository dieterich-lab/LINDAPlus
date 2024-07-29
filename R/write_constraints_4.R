write_constraints_4 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
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
  
  # IC Part
  for (cell in cell_types) {
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
      
      for(jj in seq_along(idx)) {
        ind <- which(curr_ints == interactions[jj])
        if(length(ind) > 0){
          cc1 <- c(cc1, paste0(var[idx[jj]], " - ", var[ind], " >= 0"))
          cc2 <- c(cc2, paste0(var[idx[jj]], " - ", paste0(var[ind], collapse = " - "), " <= 0"))
        }
      }
      
      constraints <- c(constraints, cc1, cc2)
    }
  }
  
  var2rem <- variables$var[which(grepl(pattern = "Cardiomyocytes:reaction PSEUDODOMAINLR=", x = variables$var_exp, fixed = TRUE))]
  if(length(var2rem) > 0){
    remrem <- c()
    for(ll in 1:length(var2rem)){
      remrem <- c(remrem, which(grepl(pattern = var2rem[ll], x = constraints, fixed = TRUE)))
    }
    remrem <- unique(remrem)
    if(length(remrem) > 0){
      constraints <- constraints[-remrem]
    }
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
  for(jj in seq_along(idx)){
    
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
