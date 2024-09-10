write_constraints_2 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  constraints <- c()
  cell_types <- names(background.networks.list$background.networks)
  
  for (cell in cell_types) {
    # Filtering and cleaning variable expressions for the current cell type
    cell_pattern <- paste0(cell, ":")
    cell_mask <- grepl(cell_pattern, variables$var_exp, fixed = TRUE)
    var <- variables$var[cell_mask]
    var_exp_raw <- variables$var_exp[cell_mask]
    var_exp <- sub(cell_pattern, "", var_exp_raw, fixed = TRUE)
    
    # Constraints for interactions
    interaction_mask <- grepl("interaction ", var_exp, fixed = TRUE)
    interaction_vars <- var[interaction_mask]
    interaction_exps <- var_exp[interaction_mask]
    lr2rem <- background.networks.list$background.networks[[cell]]
    lr2rem <- paste0("interaction ", 
                     lr2rem$gene_source[which(lr2rem$pfam_source=="PSEUDODOMAINLR")], 
                     "=", 
                     lr2rem$gene_target[which(lr2rem$pfam_source=="PSEUDODOMAINLR")])
    ind2rem <- which(interaction_exps%in%lr2rem)
    interaction_exps <- interaction_exps[-ind2rem]
    interaction_vars <- interaction_vars[-ind2rem]
    interactions <- sub("interaction ", "", interaction_exps, fixed = TRUE)
    interactions_split <- strsplit(interactions, "=")
    ss <- sapply(interactions_split, '[', 1)
    tt <- sapply(interactions_split, '[', 2)
    
    node_ss <- paste0("node ", ss)
    node_tt <- paste0("node ", tt)
    # ind2rem <- which(tt%in%background.networks.list$ligands.receptors$ligands)
    # if(length(ind2rem) > 0){
    #   node_ss <- node_ss[-ind2rem]
    #   node_tt <- node_tt[-ind2rem]
    # }
    node_ss_idx <- match(node_ss, var_exp)
    node_tt_idx <- match(node_tt, var_exp)
    
    cc1 <- paste0(var[node_ss_idx], " + ", var[node_tt_idx], " - 2 ", interaction_vars, " >= 0")
    if (length(cc1) > 0) constraints <- c(constraints, cc1)
    
    # Constraints for reactions
    reaction_mask <- grepl("reaction ", var_exp, fixed = TRUE)
    reaction_vars <- var[reaction_mask]
    reaction_exps <- var_exp[reaction_mask]
    lr2rem <- background.networks.list$background.networks[[cell]]
    lr2rem <- paste0(lr2rem$gene_source[which(lr2rem$pfam_source=="PSEUDODOMAINLR")], 
                     "=", 
                     lr2rem$gene_target[which(lr2rem$pfam_source=="PSEUDODOMAINLR")])
    lr2rem <- unique(lr2rem)
    ind2rem <- c()
    for(zz in 1:length(lr2rem)){
      ind2rem <- c(ind2rem, which(grepl(pattern = lr2rem[zz], x = reaction_exps, fixed = TRUE)))
    }
    ind2rem <- unique(ind2rem)
    reaction_exps <- reaction_exps[-ind2rem]
    reaction_vars <- reaction_vars[-ind2rem]
    reaction_details <- strsplit(reaction_exps, " ")
    ddi <- sapply(reaction_details, '[', 2)
    ppi <- sapply(reaction_details, '[', 4)
    
    dd1 <- sapply(strsplit(ddi, "="), '[', 1)
    dd2 <- sapply(strsplit(ddi, "="), '[', 2)
    pp1 <- sapply(strsplit(ppi, "="), '[', 1)
    pp2 <- sapply(strsplit(ppi, "="), '[', 2)
    
    domain1 <- paste0("domain ", dd1, " of protein ", pp1)
    domain2 <- paste0("domain ", dd2, " of protein ", pp2)
    
    domain1_idx <- match(domain1, var_exp)
    domain2_idx <- match(domain2, var_exp)
    
    cc2 <- paste0(var[domain1_idx], " + ", var[domain2_idx], " - 2 ", reaction_vars, " >= 0")
    if (length(cc2) > 0) constraints <- c(constraints, cc2)
  }
  
  return(constraints)
}
