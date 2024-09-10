write_constraints_3b <- function(variables = variables,
                                 background.networks.list = background.networks.list){
  
  constraints <- c()
  cell_types <- names(background.networks.list$background.networks)
  
  for (cell in cell_types) {
    # Filter and clean variable names for the current cell type
    cell_vars <- grepl(paste0(cell, ":"), variables$var_exp, fixed = TRUE)
    var <- variables$var[cell_vars]
    var_exp <- gsub(paste0(cell, ":"), "", variables$var_exp[cell_vars], fixed = TRUE)
    
    # Out-going constraints
    node_mask <- grepl("node ", var_exp, fixed = TRUE)
    node_vars <- var[node_mask]
    node_exp <- var_exp[node_mask]
    nodes_clean <- gsub("node ", "", node_exp, fixed = TRUE)
    
    # For each node, find all interactions where this node is the source
    cc1 <- sapply(nodes_clean, function(node) {
      interactions_mask <- grepl(paste0("interaction ", node, "="), var_exp, fixed = TRUE)
      interaction_vars <- var[interactions_mask]
      if (length(interaction_vars) > 0) {
        return(paste0(paste0(interaction_vars, collapse = " + "), " - ", var[node_mask][node == nodes_clean], " >= 0"))
      }
    })
    
    # Incoming constraints
    all_genes <- unique(c(background.networks.list$background.networks[[cell]]$gene_source,
                          background.networks.list$background.networks[[cell]]$gene_target))
    receptors <- background.networks.list$ligands.receptors$receptors
    non_receptor_nodes <- setdiff(all_genes, receptors)
    node_relevant_mask <- var_exp %in% paste0("node ", non_receptor_nodes)
    
    node_relevant_vars <- var[node_relevant_mask]
    node_relevant_exp <- var_exp[node_relevant_mask]
    node_relevant_clean <- gsub("node ", "", node_relevant_exp, fixed = TRUE)
    
    cc2 <- sapply(node_relevant_clean, function(node) {
      interactions_mask <- grepl("interaction ", var_exp, fixed = TRUE)
      interactions_second_part <- sapply(strsplit(var_exp[interactions_mask], "="), '[', 2)
      interaction_vars <- var[interactions_mask][interactions_second_part == node]
      if (length(interaction_vars) > 0) {
        return(paste0(paste0(interaction_vars, collapse = " + "), " - ", var[node_relevant_mask][node == node_relevant_clean], " >= 0"))
      }
    })
    
    cc1 <- unique(unlist(cc1))
    cc2 <- unique(unlist(cc2))
    
    tmp <- cc1[sapply(cc1, function(x) length(unlist(strsplit(x, " "))) == 5)]
    ll_exp <- paste0(cell, ":node ", background.networks.list$ligands.receptors$ligands)
    ll_var <- variables$var[which(variables$var_exp%in%ll_exp)]
    vec2rem <- vector("list", length(tmp))
    count <- 0
    for (ii in seq_along(tmp)) {
      str2check <- strsplit(tmp[ii], " ")[[1]][1]
      if (str2check %in% ll_var) {
        count <- count + 1
        vec2rem[[count]] <- tmp[ii]
      }
    }
    vec2rem <- unlist(vec2rem[1:count])

    ind2rem <- which(cc1%in%vec2rem)
    if(length(ind2rem) > 0){
      cc1 <- cc1[-ind2rem]
    }

    tmp <- cc2[sapply(cc2, function(x) length(unlist(strsplit(x, " "))) == 5)]
    ll_exp <- paste0(cell, ":node ", background.networks.list$ligands.receptors$ligands)
    ll_var <- variables$var[which(variables$var_exp%in%ll_exp)]
    vec2rem <- vector("list", length(tmp))
    count <- 0
    for (ii in seq_along(tmp)) {
      str2check <- strsplit(tmp[ii], " ")[[1]][3]
      if (str2check %in% ll_var) {
        count <- count + 1
        vec2rem[[count]] <- tmp[ii]
      }
    }
    vec2rem <- unlist(vec2rem[1:count])

    ind2rem <- which(cc2%in%vec2rem)
    if(length(ind2rem) > 0){
      cc2 <- cc2[-ind2rem]
    }
    
    # Combine and append constraints
    constraints <- c(constraints, cc1, cc2)
  }
  
  return(constraints)
  
}