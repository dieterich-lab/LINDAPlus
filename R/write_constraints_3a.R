write_constraints_3a <- function(variables = variables,
                                 background.networks.list = background.networks.list,
                                 constraits.parallel.writing = constraits.parallel.writing) {
  
  constraints <- c()
  cell_types <- names(background.networks.list$background.networks)
  
  process_cell_type <- function(cell) {
    local_constraints <- c()
    
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
      } else {
        return(NULL)
      }
    }, USE.NAMES = FALSE)
    
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
      } else {
        return(NULL)
      }
    }, USE.NAMES = FALSE)
    
    # Combine constraints
    cc1 <- unique(unlist(cc1))
    cc2 <- unique(unlist(cc2))
    
    local_constraints <- c(local_constraints, cc1, cc2)
    return(local_constraints)
  }
  
  if (constraits.parallel.writing) {
    constraints_list <- mclapply(cell_types, process_cell_type, mc.cores = length(cell_types))
  } else {
    constraints_list <- lapply(cell_types, process_cell_type)
  }
  
  constraints <- unique(unlist(constraints_list))
  
  return(constraints)
}