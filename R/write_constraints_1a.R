write_constraints_1a <- function(variables = variables,
                                 background.networks.list = background.networks.list){
  
  constraints <- c()
  
  receptors <- background.networks.list$ligands.receptors$receptors
  ligands <- background.networks.list$ligands.receptors$ligands
  
  background.network.list <- background.networks.list$background.networks
  
  # There should be at least one Receptor activated on each cell
  cell_types <- names(background.networks.list$background.networks)
  
  for(ii in 1:length(cell_types)){
    
    idx2keep <- which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])
    var_exp <- variables$var_exp[idx2keep]
    var <- variables$var[idx2keep]
    
    idx <- which(var_exp%in%paste0(cell_types[ii], ":node ", receptors))
    
    cc1 <- paste0(paste0(var[idx], collapse = " + "), " >= 1")
    
    constraints <- c(constraints, cc1)
    
  }
  
  # There should be at least one Ligand activated on each cell
  for(ii in 1:length(cell_types)){
    
    idx2keep <- which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])
    var_exp <- variables$var_exp[idx2keep]
    var <- variables$var[idx2keep]
    
    bn <- background.network.list[[cell_types[ii]]]
    
    ll <- unique(bn$gene_source[intersect(x = which(bn$pfam_source=="PSEUDODOMAINLR"), 
                                          y = which(bn$pfam_target!="PSEUDODOMAINLR"))])
    
    idx <- which(var_exp%in%paste0(cell_types[ii], ":node ", ll))
    
    if(length(idx) > 0){
      
      cc12 <- paste0(paste0(var[idx], collapse = " + "), " >= 1")
      constraints <- c(constraints, cc12)
      
    }
    
  }
  
  
  # # A receptor in cell-type k is present in the solution, then
  # # there should be at least one functional interaction with a
  # # ligand upstream of it
  # 
  # # Protein Level
  # for(ii in 1:length(cell_types)){
  #   
  #   var <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
  #   var_exp <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
  #   
  #   bg <- background.networks.list$background.networks[[ii]]
  #   
  #   lr <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR")]
  #   lr_exp <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR")]
  #   for(jj in 1:length(receptors)){
  #     
  #     idx <- intersect(x = which(bg$gene_source%in%ligands), y = which(bg$gene_target%in%receptors[jj]))
  #     if(length(idx) > 0){
  #       
  #       tmp <- bg[idx, ]
  #       ind1 <- which(var_exp==paste0(cell_types[ii], ":node ", receptors[jj]))
  #       
  #       ind2 <- which(lr_exp%in%paste0("LR:interaction ", tmp$gene_source, "=", tmp$gene_target))
  #       
  #       cc2 <- paste0(lr[ind2], collapse = " - ")
  #       cc2 <- paste0(var[ind1], " - ", cc2, " <= 0")
  #       constraints <- c(constraints, cc2)
  #       
  #       cc22 <- paste0(var[ind1], " - ", lr[ind2], " >= 0")
  #       constraints <- c(constraints, cc22)
  #       
  #     }
  #     
  #   }
  #   
  # }
  # 
  # # Domain Level
  # for(ii in 1:length(cell_types)){
  #   
  #   var <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
  #   var_exp <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
  #   
  #   bg <- background.networks.list$background.networks[[ii]]
  #   
  #   lr <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR")]
  #   lr_exp <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR")]
  #   for(jj in 1:length(receptors)){
  #     
  #     idx <- intersect(x = which(bg$gene_source%in%ligands), y = which(bg$gene_target%in%receptors[jj]))
  #     if(length(idx) > 0){
  #       
  #       tmp <- bg[idx, ]
  #       udomains <- unique(paste0(tmp$pfam_target))
  #       
  #       for(kk in 1:length(udomains)){
  #         
  #         tmptmp <- tmp[which(tmp$pfam_target==udomains[kk]), ]
  #         
  #         ind1 <- which(lr_exp%in%paste0("LR:reaction ", tmptmp$pfam_source, "_", tmptmp$gene_source,
  #                                        "=", tmptmp$pfam_target, "_", tmptmp$gene_target))
  #         
  #         ind2 <- which(var_exp==paste0(cell_types[ii], ":domain ", udomains[kk], 
  #                                       " of protein ", receptors[jj]))
  #         
  #         
  #         cc3 <- paste0(lr[ind1], collapse = " - ")
  #         cc3 <- paste0(var[ind2], " - ", cc3, " <= 0")
  #         constraints <- c(constraints, cc3)
  #         
  #         cc33 <- paste0(var[ind2], " - ", lr[ind1], " >= 0")
  #         constraints <- c(constraints, cc33)
  #         
  #       }
  #       
  #     }
  #     
  #   }
  #   
  # }
  
  
  # If there are joint receptors, then all of the proteins which they consist
  # of should be activated
  cc4 <- c()
  cc5 <- c()
  for(ii in 1:length(background.network.list)){
    
    idx <- which(grepl(pattern = "|", x = background.network.list[[ii]]$gene_source, fixed = TRUE))
    if(length(idx) > 0){
      
      complex_receptor <- unique(background.network.list[[ii]]$gene_source[idx])
      for(jj in 1:length(complex_receptor)){
        
        proteins <- unique(unlist(strsplit(x = complex_receptor[jj], split = "|", fixed = TRUE)))
        ind1 <- which(variables$var_exp==paste0(names(background.network.list)[ii], ":node ", complex_receptor[jj]))
        ind2 <- which(variables$var_exp%in%paste0(names(background.network.list)[ii], ":node ", proteins))
        
        tmp1 <- paste0(paste0(variables$var[ind2], collapse = " + "), " - ", length(ind2)-1, " ", variables$var[ind1], " >= 0")
        tmp2 <- paste0(paste0(variables$var[ind2], collapse = " + "), " - ", length(ind2)-1, " ", variables$var[ind1], " <= ", length(ind2)-1)
        
        cc4 <- c(cc4, tmp1)
        cc5 <- c(cc5, tmp2)
        
      }
      
    }
    
  }
  if((length(cc4) > 0) && (length(cc5) > 0)){
    constraints <- c(constraints, unique(cc4), unique(cc5))
  }
  
  # If a non-extra-cellular domain of a protein receptor is activated, then also
  # at least one of the extra-cellular domains of the protein receptors should
  # be activated
  cc6 <- c()
  for(ii in 1:length(cell_types)){
    
    bn <- background.networks.list$background.networks[[cell_types[ii]]]
    for(jj in 1:length(receptors)){
      
      curr <- receptors[jj]
      ind1 <- intersect(x = which(bn$gene_target==curr), 
                        y = intersect(x = which(bn$pfam_source=="PSEUDODOMAINLR"), 
                                      y = which(!grepl(pattern = "PSEUDODOMAIN", x = bn$pfam_target, fixed = TRUE))))
      ind2 <- intersect(x = which(bn$gene_target==curr), 
                        y = intersect(x = which(!grepl(pattern = "PSEUDODOMAIN", x = bn$pfam_source, fixed = TRUE)), 
                                      y = which(!grepl(pattern = "PSEUDODOMAIN", x = bn$pfam_target, fixed = TRUE))))
      
      outside_domain <- unique(bn$pfam_target[ind1])
      inside_domain <- setdiff(x = bn$pfam_target[ind2], y = outside_domain)
      
      if((length(outside_domain) > 0) && (length(inside_domain) > 0)){
        
        for(kk in 1:length(inside_domain)){
          
          ind <- intersect(x = intersect(x = which(bn$gene_target==curr), 
                                         y = which(bn$pfam_target==inside_domain[kk])), 
                           y = intersect(x = which(!grepl(pattern = "PSEUDODOMAIN", x = bn$pfam_source, fixed = TRUE)), 
                                         y = which(bn$pfam_source%in%outside_domain)))
          if(length(ind) > 0){
            
            tmp <- bn[ind, ]
            
            vv1 <- variables$var[which(variables$var_exp==paste0(cell_types[ii], ":domain ", inside_domain[kk], " of protein ", curr))]
            vv2 <- c()
            for(ll in 1:nrow(tmp)){
              vv2 <- c(vv2, variables$var[which(variables$var_exp==paste0(cell_types[ii], ":domain ", tmp$pfam_source[ll], 
                                                                          " of protein ", tmp$gene_source[ll]))])
            }
            
            cc61 <- paste0(vv1, " - ", paste0(vv2, collapse = " - "), " <= 0")
            cc62 <- paste0(vv1, " - ", vv2, " >= 0")
            cc6 <- c(cc6, c(cc61, cc62))
            
          }
          
        }
        
      }
      
    }
    
  }
  
  constraints <- unique(c(constraints, cc6))
  
  return(constraints)
  
}