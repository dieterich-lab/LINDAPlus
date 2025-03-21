process_network_solutions <- function(res = res){
  
  net <- res$combined_solutions
  cells <- unique(net$Space)
  cells <- cells[-which(cells == "Extra-Cellular")]
  lr_int <- paste0(net$Gene_Source[which(net$Space=="Extra-Cellular")], 
                   "=", 
                   net$Gene_Target[which(net$Space=="Extra-Cellular")])
  ligands <- sapply(strsplit(x = lr_int, split = "=", fixed = TRUE), "[", 1)
  receptors <- sapply(strsplit(x = lr_int, split = "=", fixed = TRUE), "[", 2)
  
  common_ligands <- intersect(x = intersect(x = net$Gene_Source[which(net$Space == "Extra-Cellular")], 
                                            y = net$Gene_Target[which(net$DDI == "PSEUDODOMAINTF=PSEUDODOMAINTF")]), 
                              y = net$Gene_Source[-which(net$Space == "Extra-Cellular")])
  
  diff_ligands <- setdiff(x = intersect(x = net$Gene_Source[which(net$Space == "Extra-Cellular")], 
                                        y = net$Gene_Target[which(net$DDI == "PSEUDODOMAINTF=PSEUDODOMAINTF")]), 
                          y = common_ligands)
  
  if(length(diff_ligands) > 0){
    
    lr_diff <- net[intersect(x = which(net$Space == "Extra-Cellular"), 
                             y = which(net$Gene_Source %in% diff_ligands)), ]
    
    tobind <- matrix(data = , nrow = 1, ncol = ncol(net))
    for(ii in 1:nrow(lr_diff)){
      
      curr_rec <- lr_diff$Gene_Target[ii]
      ind <- which(net$Gene_Source == curr_rec)
      if(length(ind) > 0){
        
        spaces <- unique(net$Space[ind])
        tobind2 <- matrix(data = , nrow = length(spaces), ncol = ncol(net))
        tobind2[, 1] <- spaces
        tobind2[, 2] <- lr_diff$Gene_Source[ii]
        tobind2[, 3] <- lr_diff$Gene_Target[ii]
        tobind2[, 4] <- lr_diff$Weight[ii]
        tobind2[, 5] <- lr_diff$DDI[ii]
        
        tobind <- unique(rbind(tobind, tobind2))
        
        
      }
      
      
    }
    tobind <- tobind[2:nrow(tobind), ]
    colnames(tobind) <- colnames(net)
    tobind <- as.data.frame(tobind)
    tobind$Weight <- as.numeric(tobind$Weight)
    net <- rbind(net, tobind)
    
    ligands <- unique(net$Gene_Source[which(net$Space == "Extra-Cellular")])
    receptors <- unique(net$Gene_Target[which(net$Space == "Extra-Cellular")])
    
    ind2rem <- unique(c(intersect(x = which(net$Gene_Source %in% ligands), y = which(net$Gene_Target %in% ligands)),
                        intersect(x = which(net$Gene_Source %in% receptors), y = which(net$Gene_Target %in% receptors)),
                        intersect(x = which(net$Gene_Source %in% receptors), y = which(net$Gene_Target %in% ligands))))
    
    if(length(ind2rem) > 0){
      net <- net[-ind2rem, ]
    }
    
    receptors <- unique(net$Gene_Target[which(net$Space == "Extra-Cellular")])
    ligands <- unique(net$Gene_Source[which(net$Space == "Extra-Cellular")])
    mm <- matrix(data = , nrow = 1, ncol = 5)
    for(ii in 1:length(receptors)){
      
      for(jj in 1:length(cells)){
        
        curr <- net[which(net$Space == cells[jj]), ]
        idx_check <- intersect(x = which(curr$Gene_Source %in% ligands), 
                               y = which(curr$Gene_Target == receptors[ii]))
        if(length(idx_check) == 0){
          
          indind <- intersect(x = which(net$Space == "Extra-Cellular"), 
                              y = which(net$Gene_Target == receptors[ii]))
          if(length(indind) > 0){
            
            tobind <- matrix(data = , nrow = length(indind), ncol = 5)
            tobind[, 1] <- cells[jj]
            tobind[, 2] <- net$Gene_Source[indind]
            tobind[, 3] <- receptors[ii]
            tobind[, 4] <- net$Weight[indind]
            tobind[, 5] <- net$DDI[indind]
            mm <- rbind(mm, tobind)
            
          }
          
        }
        
      }
      
    }
    
    if(nrow(mm) > 1){
      mm <- mm[2:nrow(mm), ]
      colnames(mm) <- colnames(net)
      mm <- as.data.frame(mm)
      mm$Weight <- as.numeric(mm$Weight)
      net <- unique(rbind(net, mm))
    }
    
    receptors <- unique(net$Gene_Target[which(net$Space == "Extra-Cellular")])
    int2rem <- c()
    for(ii in 1:length(receptors)){
      
      for(jj in 1:length(cells)){
        
        curr <- net[which(net$Space == cells[jj]), ]
        
        indind <- which(curr$Gene_Source == receptors[ii])
        if(length(indind) == 0){
          ind2rem <- which(curr$Gene_Target == receptors[ii])
          int2rem <- c(int2rem, paste0(curr$Space[ind2rem], "=", curr$Gene_Source[ind2rem], 
                                       "=", curr$Gene_Target[ind2rem]))
        }
        
      }
      
    }
    if(length(int2rem) > 0){
      ind2rem <- which(paste0(net$Space, "=", net$Gene_Source, "=", net$Gene_Target) %in% int2rem)
      if(length(ind2rem) > 0){
        net <- net[-ind2rem, ]
      }
    }
    
    net_ec <- net[which(net$Space == "Extra-Cellular"), ]
    net_other <- net[which(net$Space != "Extra-Cellular"), ]
    lr2rem <- setdiff(x = unique(paste0(net_ec$Gene_Source, "=", net_ec$Gene_Target)), 
                      y = unique(paste0(net_other$Gene_Source, "=", net_other$Gene_Target)))
    if(length(lr2rem) > 0){
      ind2rem <- which(paste0(net$Space, "=", net$Gene_Source, "=", net$Gene_Target) %in% 
                         paste0("Extra-Cellular=", lr2rem))
      if(length(ind2rem) > 0){
        net <- net[-ind2rem, ]
      }
    }
    
    
    #### Process Attributes table
    atr <- res$node_attributes
    
    all_network_genes <- unique(c(net$Gene_Source, net$Gene_Target))
    all_attributes_genes <- sapply(strsplit(x = atr$node, split = "_", fixed = TRUE), "[", 1)
    
    ind2keep <- which(all_attributes_genes %in% all_network_genes)
    if(length(ind2keep) > 0){
      
      atr <- atr[ind2keep, ]
      
    }
    
    #### Process free receptors
    extract_pairs <- function(path) {
      genes <- names(path)
      if (length(genes) > 1) {
        from <- genes[-length(genes)]   # All but the last gene
        to <- genes[-1]                 # All but the first gene
        data.frame(From = from, To = to)
      } else {
        NULL  # In case there's a path with a single vertex, return nothing
      }
    }
    free_receptors <- setdiff(x = atr$node[which(atr$attribute == "Receptor")], 
                              y = net$Gene_Target[which(net$Space == "Extra-Cellular")])
    if(length(free_receptors) > 0){
      
      ind2rem <- c()
      for(ii in 1:length(free_receptors)){
        
        cells <- unique(net$Space[which(net$Gene_Source == free_receptors[ii])])
        for(jj in 1:length(cells)){
          
          curr <- net[which(net$Space == cells[jj]), ]
          df <- unique(curr[, 2:3])
          gg <- igraph::graph_from_data_frame(d = df, directed = TRUE)
          adj <- igraph::as_adjacency_matrix(graph = gg)
          tf <- intersect(x = colnames(adj), y = atr$node[which(atr$attribute == "TF")])
          receptors <- setdiff(x = intersect(x = colnames(adj), y = net$Gene_Target[which(net$Space == "Extra-Cellular")]), 
                               y = tf)
          if((length(tf) > 0) && (length(receptors) > 0)){
            
            # sp_free <- igraph::all_simple_paths(graph = gg, from = which(rownames(adj) == free_receptors[ii]), to = which(rownames(adj) %in% tf))
            sp_true <- list()
            for(kk in 1:length(receptors)){
              sp <- igraph::all_simple_paths(graph = gg, from = which(rownames(adj) == receptors[kk]), to = which(rownames(adj) %in% tf))
              for(ll in 1:length(sp)){
                sp_true[[length(sp_true)+1]] <- sp[[ll]]
              }
            }
            
            # df_free <- do.call(rbind, lapply(sp_free, extract_pairs))
            df_true <- do.call(rbind, lapply(sp_true, extract_pairs))
            int2rem <- setdiff(x = paste0(curr$Gene_Source, "=", curr$Gene_Target), 
                               y = paste0(df_true$From, "=", df_true$To))
            indind <- which(sapply(strsplit(x = int2rem, split = "=", fixed = TRUE), "[", 2) %in% receptors)
            if(length(indind) > 0){int2rem <- int2rem[-indind]}
            if(length(int2rem) > 0){
              ind2rem <- c(ind2rem, which(paste0(net$Space, "=", net$Gene_Source, "=", net$Gene_Target) %in%
                                            paste0(cells[jj], "=", int2rem)))
            }
            
          }
          
        }
        
        
      }
      
      if(length(ind2rem) > 0){
        
        net <- net[-ind2rem, ]
        
      }
      
    }
    
    all_network_genes <- unique(c(net$Gene_Source, net$Gene_Target))
    all_attributes_genes <- sapply(strsplit(x = atr$node, split = "_", fixed = TRUE), "[", 1)
    
    ind2keep <- which(all_attributes_genes %in% all_network_genes)
    if(length(ind2keep) > 0){
      
      atr <- atr[ind2keep, ]
      
    }
    
    res <- list()
    res[[length(res)+1]] <- net
    res[[length(res)+1]] <- atr
    names(res) <- c("combined_solutions", "node_attributes")
    
    return(res)
    
  } else {
    
    return(res)
    
  }
  
}