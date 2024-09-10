create_variables <- function(background.networks.list = background.networks.list){
  
  background.network.list <- background.networks.list$background.networks
  
  cc1 <- c()
  cc1Exp <- c()
  for(ii in 1:length(background.network.list)){
    
    background.network <- background.network.list[[ii]]
    
    cc1 <- c(cc1, paste0("xb", 1:nrow(background.network), "_", names(background.network.list)[ii]))
    cc1Exp <- c(cc1Exp, paste0(names(background.network.list)[ii], ":",
                               "reaction ", background.network$pfam_source, "=",
                               background.network$pfam_target, " of ",
                               background.network$gene_source, "=",
                               background.network$gene_target))
  }
  
  
  #
  cc2 <- c()
  cc2Exp <- c()
  for(ii in 1:length(background.network.list)){
    
    cnt <- 1 + length(which(sapply(strsplit(x = cc1, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii]))
    ppi <- unique(background.network.list[[ii]][, c("gene_source", "gene_target")])
    cc2 <- c(cc2, paste0("xb", cnt:(nrow(ppi)+cnt-1), "_", names(background.network.list)[ii]))
    cc2Exp <- c(cc2Exp, paste0(names(background.network.list)[ii], ":interaction ", 
                               ppi$gene_source, "=", ppi$gene_target))
    
  }
  
  
  #
  cc3 <- c()
  cc3Exp <- c()
  for(ii in 1:length(background.network.list)){
    
    cnt <- 1 + length(which(sapply(strsplit(x = cc1, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii])) +
      length(which(sapply(strsplit(x = cc2, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii]))
    ppi <- unique(background.network.list[[ii]][, c("gene_source", "gene_target")])
    nodes <- unique(c(ppi$gene_source, ppi$gene_target))
    cc3 <- c(cc3, paste0("xb", cnt:(length(nodes)+cnt-1), "_", names(background.network.list)[ii]))
    cc3Exp <- c(cc3Exp, paste0(names(background.network.list)[ii], ":node ", nodes))
    
  }
  
  
  #
  cc4 <- c()
  cc4Exp <- c()
  for(ii in 1:length(background.network.list)){
    
    cnt <- 1 + length(which(sapply(strsplit(x = cc1, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii])) +
      length(which(sapply(strsplit(x = cc2, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii])) +
      length(which(sapply(strsplit(x = cc3, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii]))
    ppi <- unique(background.network.list[[ii]])
    tmp <- unique(c(paste0(names(background.network.list)[ii], 
                           ":domain ",
                           ppi$pfam_source,
                           " of protein ",
                           ppi$gene_source),
                    paste0(names(background.network.list)[ii],
                           ":domain ",
                           ppi$pfam_target,
                           " of protein ",
                           ppi$gene_target)))
    cc4Exp <- c(cc4Exp, tmp)
    cc4 <- c(cc4, paste0("xb", cnt:(length(tmp)+cnt-1), 
                         "_", names(background.network.list)[ii]))
    
  }
  
  
  #
  cc5 <- c()
  cc5Exp <- c()
  for(ii in 1:length(background.network.list)){
    
    cnt <- 1 + length(which(sapply(strsplit(x = cc1, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii])) +
      length(which(sapply(strsplit(x = cc2, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii])) +
      length(which(sapply(strsplit(x = cc3, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii])) +
      length(which(sapply(strsplit(x = cc4, split = "_", fixed = TRUE), "[", 2)==names(background.network.list)[ii]))
    ppi <- unique(background.network.list[[ii]][, c("gene_source", "gene_target")])
    nodes <- unique(c(ppi$gene_source, ppi$gene_target))
    cc5 <- c(cc5, paste0("xb", cnt:(length(nodes)+cnt-1), "_", names(background.network.list)[ii]))
    cc5Exp <- c(cc5Exp, paste0(names(background.network.list)[ii], ":dist ", nodes))
    
  }
  
  
  #
  cnt <- 1
  for(ii in 1:length(background.network.list)){
    
    bg <- background.network.list[[ii]]
    
    idx2keep <- intersect(x = which(bg$gene_source%in%background.networks.list$ligands.receptors$ligand), 
                          y = which(bg$gene_target%in%background.networks.list$ligands.receptors$receptor))
    
    if((length(idx2keep) > 0) && (cnt == 1)){
      lr <- bg[idx2keep, ]
      cnt <- cnt + 1
    } else {
      if(length(idx2keep) > 0){
        lr <- unique(rbind(lr, bg[idx2keep, ]))
      }
    }
  }
  
  receptors <- paste0("receptor ", unique(intersect(x = lr$gene_target, y = background.networks.list$ligands.receptors$receptors)))
  ligands <- paste0("ligand ", unique(intersect(x = lr$gene_source, y = background.networks.list$ligands.receptors$ligands)))
  
  receptors_domains <- paste0("domain_receptor ", unique(paste0(lr$pfam_target, "_", lr$gene_target)))
  ligands_domains <- paste0("domain_ligand ", unique(paste0(lr$pfam_source, "_", lr$gene_source)))
  
  reactions <- paste0("reaction ", unique(paste0(lr$pfam_source, "_", lr$gene_source, "=", lr$pfam_target, "_", lr$gene_target)))
  interactions <- paste0("interaction ", unique(paste0(lr$gene_source, "=", lr$gene_target)))
  
  cc6Exp <- paste0("LR:", c(receptors, ligands, receptors_domains, ligands_domains, reactions, interactions))
  cc6 <- paste0("lr", 1:length(cc6Exp))
  
  
  #
  cc <- c(cc1, cc2, cc3, cc4, cc5, cc6)
  ccExp <- c(cc1Exp, cc2Exp, cc3Exp, cc4Exp, cc5Exp, cc6Exp)
  
  variables <- list()
  variables[[length(variables)+1]] <- cc
  variables[[length(variables)+1]] <- ccExp
  names(variables) <- c("var", "var_exp")
  return(variables)
  
}
