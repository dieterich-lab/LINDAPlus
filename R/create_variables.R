create_variables <- function(background.networks.list = background.networks.list,
                             constraits.parallel.writing = constraits.parallel.writing){
  
  background.network.list <- background.networks.list$background.networks
  num_cores <- length(background.network.list)
  
  process_network <- function(ii) {
    background.network <- background.network.list[[ii]]
    network_name <- names(background.network.list)[ii]
    
    # cc1 and cc1Exp
    cc1 <- paste0("xb", 1:nrow(background.network), "_", network_name)
    cc1Exp <- paste0(network_name, ":reaction ", background.network$pfam_source, "=", background.network$pfam_target,
                     " of ", background.network$gene_source, "=", background.network$gene_target)
    
    # cc2 and cc2Exp
    cnt <- 1 + sum(sapply(strsplit(cc1, "_", fixed = TRUE), "[", 2) == network_name)
    ppi <- unique(background.network[, c("gene_source", "gene_target")])
    cc2 <- paste0("xb", cnt:(nrow(ppi) + cnt - 1), "_", network_name)
    cc2Exp <- paste0(network_name, ":interaction ", ppi$gene_source, "=", ppi$gene_target)
    
    # cc3 and cc3Exp
    cnt <- cnt + length(cc2)
    nodes <- unique(c(ppi$gene_source, ppi$gene_target))
    cc3 <- paste0("xb", cnt:(length(nodes) + cnt - 1), "_", network_name)
    cc3Exp <- paste0(network_name, ":node ", nodes)
    
    # cc4 and cc4Exp
    cnt <- cnt + length(cc3)
    tmp <- unique(c(paste0(network_name, ":domain ", background.network$pfam_source, " of protein ", background.network$gene_source),
                    paste0(network_name, ":domain ", background.network$pfam_target, " of protein ", background.network$gene_target)))
    cc4 <- paste0("xb", cnt:(length(tmp) + cnt - 1), "_", network_name)
    cc4Exp <- tmp
    
    # cc5 and cc5Exp
    cnt <- cnt + length(cc4)
    cc5 <- paste0("xb", cnt:(length(nodes) + cnt - 1), "_", network_name)
    cc5Exp <- paste0(network_name, ":dist ", nodes)
    
    list(cc1, cc1Exp, cc2, cc2Exp, cc3, cc3Exp, cc4, cc4Exp, cc5, cc5Exp)
  }
  
  if (constraits.parallel.writing) {
    results <- mclapply(1:num_cores, process_network, mc.cores = num_cores)
  } else {
    results <- lapply(1:num_cores, process_network)
  }
  
  # Extract and combine results
  cc1 <- unlist(lapply(results, "[[", 1))
  cc1Exp <- unlist(lapply(results, "[[", 2))
  cc2 <- unlist(lapply(results, "[[", 3))
  cc2Exp <- unlist(lapply(results, "[[", 4))
  cc3 <- unlist(lapply(results, "[[", 5))
  cc3Exp <- unlist(lapply(results, "[[", 6))
  cc4 <- unlist(lapply(results, "[[", 7))
  cc4Exp <- unlist(lapply(results, "[[", 8))
  cc5 <- unlist(lapply(results, "[[", 9))
  cc5Exp <- unlist(lapply(results, "[[", 10))
  
  # Process ligand-receptor pairs
  cnt <- 1
  lr <- NULL
  for(ii in 1:length(background.network.list)){
    bg <- background.network.list[[ii]]
    idx2keep <- intersect(which(bg$gene_source %in% background.networks.list$ligands.receptors$ligand),
                          which(bg$gene_target %in% background.networks.list$ligands.receptors$receptor))
    
    if(length(idx2keep) > 0) {
      if(cnt == 1) {
        lr <- bg[idx2keep, ]
        cnt <- cnt + 1
      } else {
        lr <- unique(rbind(lr, bg[idx2keep, ]))
      }
    }
  }
  
  receptors <- paste0("receptor ", unique(intersect(lr$gene_target, background.networks.list$ligands.receptors$receptors)))
  ligands <- paste0("ligand ", unique(intersect(lr$gene_source, background.networks.list$ligands.receptors$ligands)))
  receptors_domains <- paste0("domain_receptor ", unique(paste0(lr$pfam_target, "_", lr$gene_target)))
  ligands_domains <- paste0("domain_ligand ", unique(paste0(lr$pfam_source, "_", lr$gene_source)))
  reactions <- paste0("reaction ", unique(paste0(lr$pfam_source, "_", lr$gene_source, "=", lr$pfam_target, "_", lr$gene_target)))
  interactions <- paste0("interaction ", unique(paste0(lr$gene_source, "=", lr$gene_target)))
  
  cc6Exp <- paste0("LR:", c(receptors, ligands, receptors_domains, ligands_domains, reactions, interactions))
  cc6 <- paste0("lr", 1:length(cc6Exp))
  
  # Combine all variables
  cc <- c(cc1, cc2, cc3, cc4, cc5, cc6)
  ccExp <- c(cc1Exp, cc2Exp, cc3Exp, cc4Exp, cc5Exp, cc6Exp)
  
  variables <- list(var = cc, var_exp = ccExp)
  return(variables)
}
