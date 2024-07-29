write_constraints_9 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  constraints <- c()
  
  # Remove the LR interactions that are present as PPI/DDI in DIGGER
  for(ii in 1:length(background.networks.list)){
    
    curr_cell <- names(background.networks.list$background.networks)[ii]
    
    bn <- background.networks.list$background.networks[[ii]]
    ligands <- intersect(x = background.networks.list$ligands.receptors$ligands, 
                         y = bn$gene_source[which(bn$pfam_source=="PSEUDODOMAINLR")])
    
    reacs <- intersect(x = variables$var_exp[which(grepl(pattern = "reaction", x = variables$var_exp, fixed = TRUE))], 
                       y = variables$var_exp[which(grepl(pattern = paste0(curr_cell, ":"), x = variables$var_exp, fixed = TRUE))])
    
    ind2keep <- c()
    for(jj in 1:length(ligands)){
      idx <- which(grepl(pattern = paste0(" of ", ligands[jj], "="), x = reacs, fixed = TRUE))
      if(length(idx) > 0){ind2keep <- c(ind2keep, idx)}
    }
    reacs <- reacs[unique(ind2keep)]
    reacs <- reacs[-which(grepl(pattern = "PSEUDODOMAIN", x = reacs, fixed = TRUE))]
    
    if(length(reacs) > 0){
      constraints <- c(constraints,
                       paste0(variables$var[which(variables$var_exp%in%reacs)], " = 0"))
    }
    
  }
  
  # Now also remove any possible LR interactions involving a ligand which is not
  # being produced by any TF.
  
  ligands_lr <- c()
  ligands_tf <- c()
  for(ii in 1:length(background.networks.list$background.networks)){
    curr <- background.networks.list$background.networks[[ii]]
    ligands_lr <- c(ligands_lr, curr$gene_source[which(curr$pfam_source=="PSEUDODOMAINLR")])
    ligands_tf <- c(ligands_tf, curr$gene_target[which(curr$pfam_source=="PSEUDODOMAINTF")])
  }
  ligands_lr <- unique(ligands_lr)
  ligands_tf <- unique(ligands_tf)
  
  lig2rem <- setdiff(x = ligands_lr, y = ligands_tf)
  idxidx <- which(lig2rem == "PSEUDOLIGAND")
  if(length(idxidx) > 0){lig2rem <- lig2rem[-idxidx]}
  if(length(lig2rem) > 0){
    
    var2rem <- c()
    for(ii in 1:length(lig2rem)){
      
      curr_lig <- lig2rem[ii]
      for(jj in 1:length(background.networks.list$background.networks)){
        
        curr <- background.networks.list$background.networks[[jj]]
        
        ind <- intersect(x = which(curr$pfam_source=="PSEUDODOMAINLR"), 
                         y = which(curr$gene_source==curr_lig))
        
        if(length(ind) > 0){
          
          reacs <- paste0(names(background.networks.list$background.networks)[jj], ":",
                          "reaction ", curr$pfam_source[ind], "=", curr$pfam_target[ind],
                          " of ", curr$gene_source[ind], "=", curr$gene_target[ind])
          
          idx <- which(variables$var_exp %in% reacs)
          if(length(idx) > 0){
            var2rem <- c(var2rem, variables$var[idx])
          }
          
        }
        
      }
      
    }
    
    var2rem <- unique(var2rem)
    
    if(length(var2rem) > 0){
      constraints <- c(constraints, paste0(var2rem, " = 0"))
    }
    
  }
  
  return(constraints)
  
}