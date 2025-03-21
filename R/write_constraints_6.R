write_constraints_6 <- function(variables = variables,
                                background.networks.list = background.networks.list,
                                constraits.parallel.writing = constraits.parallel.writing) {
  
  constraints <- c()
  
  # Extract reaction and interaction variables
  exp_reactions <- variables$var_exp[grepl("LR:reaction ", variables$var_exp, fixed = TRUE)]
  var_reactions <- variables$var[grepl("LR:reaction ", variables$var_exp, fixed = TRUE)]
  
  exp_interactions <- variables$var_exp[grepl("LR:interaction ", variables$var_exp, fixed = TRUE)]
  var_interactions <- variables$var[grepl("LR:interaction ", variables$var_exp, fixed = TRUE)]
  
  # Process reactions
  process_reactions <- function(ii) {
    curr <- gsub("LR:reaction ", "", exp_reactions[ii], fixed = TRUE)
    sDomain <- strsplit(strsplit(curr, "=", fixed = TRUE)[[1]][1], "_", fixed = TRUE)[[1]][1]
    sGene <- strsplit(strsplit(curr, "=", fixed = TRUE)[[1]][1], "_", fixed = TRUE)[[1]][2]
    tDomain <- strsplit(strsplit(curr, "=", fixed = TRUE)[[1]][2], "_", fixed = TRUE)[[1]][1]
    tGene <- strsplit(strsplit(curr, "=", fixed = TRUE)[[1]][2], "_", fixed = TRUE)[[1]][2]
    
    idx <- which(grepl(paste0("reaction ", sDomain, "=", tDomain, " of ", sGene, "=", tGene), 
                       variables$var_exp, fixed = TRUE))
    
    cc1 <- paste0(paste(variables$var[idx], collapse = " + "), " - ", var_reactions[ii], " >= 0")
    cc2 <- paste0(variables$var[idx], " - ", var_reactions[ii], " <= 0")
    
    return(c(cc1, cc2))
  }
  
  # Process interactions
  process_interactions <- function(ii) {
    curr <- gsub("LR:interaction ", "", exp_interactions[ii], fixed = TRUE)
    idx <- setdiff(which(grepl(paste0("interaction ", curr), variables$var_exp, fixed = TRUE)), 
                   which(grepl("LR:interaction ", variables$var_exp, fixed = TRUE)))
    
    cc1 <- paste0(paste(variables$var[idx], collapse = " + "), " - ", var_interactions[ii], " >= 0")
    cc2 <- paste0(variables$var[idx], " - ", var_interactions[ii], " <= 0")
    
    return(c(cc1, cc2))
  }
  
  if (constraits.parallel.writing) {
    cc_reactions <- mclapply(seq_along(exp_reactions), process_reactions, mc.cores = length(background.networks.list$background.networks))
    cc_interactions <- mclapply(seq_along(exp_interactions), process_interactions, mc.cores = length(background.networks.list$background.networks))
  } else {
    cc_reactions <- lapply(seq_along(exp_reactions), process_reactions)
    cc_interactions <- lapply(seq_along(exp_interactions), process_interactions)
  }
  
  # Process zero-weight interactions
  process_zero_weight <- function(ii) {
    curr <- background.networks.list$background.networks[[ii]]
    idx <- which(curr$weight == 0)
    if (length(idx) > 0) {
      uint <- unique(paste0(curr$gene_source[idx], "=", curr$gene_target[idx]))
      ind <- which(variables$var_exp %in% paste0(names(background.networks.list$background.networks)[ii], ":interaction ", uint))
      return(paste0(variables$var[ind], " = 0"))
    }
    return(NULL)
  }
  
  if (constraits.parallel.writing) {
    cc_zero_weight <- mclapply(seq_along(background.networks.list$background.networks), process_zero_weight, mc.cores = length(background.networks.list$background.networks))
  } else {
    cc_zero_weight <- lapply(seq_along(background.networks.list$background.networks), process_zero_weight)
  }
  
  constraints <- unique(c(unlist(cc_reactions), unlist(cc_interactions), unlist(cc_zero_weight)))
  
  return(constraints)
}
