write_constraints_6 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  constraints <- c()
  
  ### Write the first part of the constrait where it creates a relationship
  ### between LR interactions in the extra-cellular space with those specified
  ### for each cell-type
  exp_reactions <- variables$var_exp[which(grepl(pattern = "LR:reaction ", x = variables$var_exp, fixed = TRUE))]
  var_reactions <- variables$var[which(grepl(pattern = "LR:reaction ", x = variables$var_exp, fixed = TRUE))]
  
  exp_interactions <- variables$var_exp[which(grepl(pattern = "LR:interaction ", x = variables$var_exp, fixed = TRUE))]
  var_interactions <- variables$var[which(grepl(pattern = "LR:interaction ", x = variables$var_exp, fixed = TRUE))]
  
  cc1 <- c()
  cc2 <- c()
  
  for(ii in 1:length(exp_reactions)){
    
    curr <- gsub(pattern = "LR:reaction ", replacement = "", x = exp_reactions[ii], fixed = TRUE)
    sDomain <- strsplit(x = strsplit(x = curr, split = "=", fixed = TRUE)[[1]][1], split = "_", fixed = TRUE)[[1]][1]
    sGene <- strsplit(x = strsplit(x = curr, split = "=", fixed = TRUE)[[1]][1], split = "_", fixed = TRUE)[[1]][2]
    tDomain <- strsplit(x = strsplit(x = curr, split = "=", fixed = TRUE)[[1]][2], split = "_", fixed = TRUE)[[1]][1]
    tGene <- strsplit(x = strsplit(x = curr, split = "=", fixed = TRUE)[[1]][2], split = "_", fixed = TRUE)[[1]][2]
    
    
    idx <- which(grepl(pattern = paste0("reaction ", sDomain, "=", tDomain, " of ", sGene, "=", tGene), 
                       x = variables$var_exp, 
                       fixed = TRUE))
    
    cc1 <- c(cc1, paste0(paste0(variables$var[idx], collapse = " + "), " - ", var_reactions[ii], " >= 0"))
    cc2 <- c(cc2, paste0(variables$var[idx], " - ", var_reactions[ii], " <= 0"))
    
  }
  
  for(ii in 1:length(exp_interactions)){
    
    curr <- gsub(pattern = "LR:interaction ", replacement = "", x = exp_interactions[ii], fixed = TRUE)
    
    idx <- setdiff(x = which(grepl(pattern = paste0("interaction ", curr), x = variables$var_exp, fixed = TRUE)), 
                   y = which(grepl(pattern = "LR:interaction ", x = variables$var_exp, fixed = TRUE)))
    
    cc1 <- c(cc1, paste0(paste0(variables$var[idx], collapse = " + "), " - ", var_interactions[ii], " >= 0"))
    cc2 <- c(cc2, paste0(variables$var[idx], " - ", var_interactions[ii], " <= 0"))
    
  }
  
  
  ### Write the constraints involving interactions with 0 probability
  cc3 <- c()
  background.network.list <- background.networks.list$background.networks
  for(ii in 1:length(background.network.list)){
    
    curr <- background.network.list[[ii]]
    idx <- which(curr$weight==0)
    if(length(idx) > 0){
      
      uint <- unique(paste0(curr$gene_source[idx], "=", curr$gene_target[idx]))
      ind <- which(variables$var_exp%in%paste0(names(background.network.list)[ii], ":interaction ", uint))
      cc3 <- c(cc3, paste0(variables$var[ind], " = 0"))
      
    }
    
  }
  
  constraints <- c(unique(cc1), unique(cc2), unique(cc3))
  
  return(constraints)
  
}