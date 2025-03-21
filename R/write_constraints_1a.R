write_constraints_1a <- function(variables = variables,
                                 background.networks.list = background.networks.list,
                                 constraits.parallel.writing = constraits.parallel.writing) {
  
  constraints <- c()
  
  receptors <- background.networks.list$ligands.receptors$receptors
  ligands <- background.networks.list$ligands.receptors$ligands
  background.network.list <- background.networks.list$background.networks
  cell_types <- names(background.networks.list$background.networks)
  
  process_cell_type <- function(ii) {
    local_constraints <- c()
    
    # There should be at least one Receptor activated on each cell
    idx2keep <- which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1) == cell_types[ii])
    var_exp <- variables$var_exp[idx2keep]
    var <- variables$var[idx2keep]
    idx <- which(var_exp %in% paste0(cell_types[ii], ":node ", receptors))
    cc1 <- paste0(paste0(var[idx], collapse = " + "), " >= 1")
    local_constraints <- c(local_constraints, cc1)
    
    # There should be at least one Ligand activated on each cell
    bn <- background.network.list[[cell_types[ii]]]
    ll <- unique(bn$gene_source[intersect(which(bn$pfam_source == "PSEUDODOMAINLR"), which(bn$pfam_target != "PSEUDODOMAINLR"))])
    idx <- which(var_exp %in% paste0(cell_types[ii], ":node ", ll))
    if (length(idx) > 0) {
      cc12 <- paste0(paste0(var[idx], collapse = " + "), " >= 1")
      local_constraints <- c(local_constraints, cc12)
    }
    
    # If a non-extra-cellular domain of a protein receptor is activated, then also at least one of the extra-cellular domains should be activated
    cc6 <- c()
    for (jj in 1:length(receptors)) {
      curr <- receptors[jj]
      ind1 <- intersect(which(bn$gene_target == curr), intersect(which(bn$pfam_source == "PSEUDODOMAINLR"), which(!grepl("PSEUDODOMAIN", bn$pfam_target))))
      ind2 <- intersect(which(bn$gene_target == curr), intersect(which(!grepl("PSEUDODOMAIN", bn$pfam_source)), which(!grepl("PSEUDODOMAIN", bn$pfam_target))))
      
      outside_domain <- unique(bn$pfam_target[ind1])
      inside_domain <- setdiff(bn$pfam_target[ind2], outside_domain)
      
      if (length(outside_domain) > 0 && length(inside_domain) > 0) {
        for (kk in 1:length(inside_domain)) {
          ind <- intersect(intersect(which(bn$gene_target == curr), which(bn$pfam_target == inside_domain[kk])), 
                           intersect(which(!grepl("PSEUDODOMAIN", bn$pfam_source)), which(bn$pfam_source %in% outside_domain)))
          if (length(ind) > 0) {
            tmp <- bn[ind, ]
            vv1 <- variables$var[which(variables$var_exp == paste0(cell_types[ii], ":domain ", inside_domain[kk], " of protein ", curr))]
            vv2 <- c()
            for (ll in 1:nrow(tmp)) {
              vv2 <- c(vv2, variables$var[which(variables$var_exp == paste0(cell_types[ii], ":domain ", tmp$pfam_source[ll], " of protein ", tmp$gene_source[ll]))])
            }
            cc61 <- paste0(vv1, " - ", paste0(vv2, collapse = " - "), " <= 0")
            cc62 <- paste0(vv1, " - ", vv2, " >= 0")
            cc6 <- c(cc6, c(cc61, cc62))
          }
        }
      }
    }
    
    local_constraints <- unique(c(local_constraints, cc6))
    return(local_constraints)
  }
  
  if (constraits.parallel.writing) {
    constraints_list <- mclapply(1:length(cell_types), process_cell_type, mc.cores = length(cell_types))
  } else {
    constraints_list <- lapply(1:length(cell_types), process_cell_type)
  }
  
  constraints <- unique(unlist(constraints_list))
  
  # Handle joint receptors
  cc4 <- c()
  cc5 <- c()
  for (ii in 1:length(background.network.list)) {
    idx <- which(grepl("|", background.network.list[[ii]]$gene_source, fixed = TRUE))
    if (length(idx) > 0) {
      complex_receptor <- unique(background.network.list[[ii]]$gene_source[idx])
      for (jj in 1:length(complex_receptor)) {
        proteins <- unique(unlist(strsplit(complex_receptor[jj], "|", fixed = TRUE)))
        ind1 <- which(variables$var_exp == paste0(names(background.network.list)[ii], ":node ", complex_receptor[jj]))
        ind2 <- which(variables$var_exp %in% paste0(names(background.network.list)[ii], ":node ", proteins))
        
        tmp1 <- paste0(paste0(variables$var[ind2], collapse = " + "), " - ", length(ind2) - 1, " ", variables$var[ind1], " >= 0")
        tmp2 <- paste0(paste0(variables$var[ind2], collapse = " + "), " - ", length(ind2) - 1, " ", variables$var[ind1], " <= ", length(ind2) - 1)
        
        cc4 <- c(cc4, tmp1)
        cc5 <- c(cc5, tmp2)
      }
    }
  }
  
  if (length(cc4) > 0 && length(cc5) > 0) {
    constraints <- unique(c(constraints, cc4, cc5))
  }
  
  return(constraints)
}
