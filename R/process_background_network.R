process_background_network <- function(background.networks.list = background.networks.list, 
                                       tf.scores = tf.scores,
                                       lr.scores = lr.scores,
                                       ccc.scores = ccc.scores){
  
  ### Prepare core BN for each cell-type
  bnList <- background.networks.list$background.networks
  temp <- list()
  for(ii in 1:length(bnList)){
    
    bn <- bnList[[ii]]
    
    idx <- intersect(x = which(is.na(bn$pfam_source)), 
                     y = which(is.na(bn$pfam_target)))
    
    ddi <- matrix(data = , nrow = 1, ncol = 4)
    if(length(idx) > 0){
      
      for(jj in 1:length(idx)){
        
        domains <- c(bn$pfam_source[which(bn$gene_source==bn$gene_source[idx[jj]])],
                     bn$pfam_target[which(bn$gene_target==bn$gene_source[idx[jj]])])
        domains <- domains[complete.cases(domains)]
        
        if(length(domains) > 0){
          
          tobind <- matrix(data = , nrow = length(domains), ncol = 4)
          tobind[, 1] <- domains
          tobind[, 2] <- "PSEUDODOMAINTF"
          tobind[, 3] <- bn$gene_source[idx[jj]]
          tobind[, 4] <- bn$gene_target[idx[jj]]
          
          ddi <- rbind(ddi, tobind)
          
        }
        
      }
      
    }
    
    if(nrow(ddi) > 1){
      if(nrow(ddi) > 2){
        ddi <- ddi[2:nrow(ddi), ]
        colnames(ddi) <- colnames(bn)
        ddi <- as.data.frame(ddi)
        bn <- rbind(bn, ddi)
      } else {
        ddi <- as.matrix(ddi[2:nrow(ddi), ])
        ddi <- t(ddi)
        colnames(ddi) <- colnames(bn)
        ddi <- as.data.frame(ddi)
        bn <- rbind(bn, ddi)
      }
    }
    
    idx2rem <- intersect(x = which(is.na(bn$pfam_source)), 
                         y = which(is.na(bn$pfam_target)))
    if(length(idx2rem) > 0){
      bn <- bn[-idx2rem, ]
    }
    
    idx <- which(is.na(bn$pfam_source))
    if(length(idx) > 0){
      bn$pfam_source[idx] <- "PSEUDODOMAINLR"
    }
    
    idx <- which(is.na(bn$pfam_target))
    if(length(idx) > 0){
      bn$pfam_target[idx] <- "PSEUDODOMAIN"
    }
    
    bn <- bn[complete.cases(bn), ]
    
    idx <- intersect(x = which(bn$gene_source%in%tf.scores[[ii]]$tf), 
                     y = which(bn$gene_target%in%background.networks.list$ligands.receptors$ligands))
    
    if(length(idx) > 0){
      bn$pfam_source[idx] <- "PSEUDODOMAINTF"
      bn$pfam_target[idx] <- "PSEUDODOMAINTF"
    }
    
    temp[[length(temp)+1]] <- bn
    
  }
  names(temp) <- names(bnList)
  
  
  ### Do LR enrichment score assignment
  temp2 <- list()
  for(ii in 1:length(temp)){
    
    df <- temp[[ii]]
    ss <- rep(NA, nrow(df))
    idx <- intersect(x = which(df$pfam_source=="PSEUDODOMAINLR"), 
                     y = which(df$pfam_target!="PSEUDODOMAINLR"))
    ss[idx] <- NA
    df$lr.scores <- ss
    
    temp2[[length(temp2)+1]] <- unique(df)
    
  }
  names(temp2) <- names(temp)
  
  if(!is.null(lr.scores)){
    
    for(ii in 1:length(lr.scores)){
      
      ind <- which(names(temp2)==names(lr.scores)[ii])
      
      if(length(ind) > 0){
        
        lr <- lr.scores[[ii]]
        # idxidx <- which(lr$score==1)
        # if(length(idxidx) > 0){
        #   lr$score[idxidx] <- 0.999
        # }
        for(jj in 1:nrow(lr)){
          
          ss <- strsplit(x = lr$lr.interaction[jj], split = "=", fixed = TRUE)[[1]][1]
          tt <- strsplit(x = lr$lr.interaction[jj], split = "=", fixed = TRUE)[[1]][2]
          
          idx <- intersect(x = which(temp2[[ind]]$gene_source==ss), 
                           y = which(temp2[[ind]]$gene_target==tt))
          if(length(idx) > 0){
            
            temp2[[ind]]$lr.scores[idx] <- abs(lr$score[jj])
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  ### Do the CCC score assignments
  if(is.null(ccc.scores)){
    
    temp3 <- vector(mode = "list", length = length(background.networks.list$background.networks))
    names(temp3) <- names(background.networks.list$background.networks)
    for(ii in 1:length(temp3)){
      temp3[[ii]] <- rep(NA, nrow(temp2[[ii]]))
    }
    
  } else {
    
    temp3 <- vector(mode = "list", length = length(background.networks.list$background.networks))
    names(temp3) <- names(background.networks.list$background.networks)
    for(ii in 1:length(temp3)){
      temp3[[ii]] <- rep(NA, nrow(temp2[[ii]]))
    }
    
    for(ii in 1:nrow(ccc.scores)){
      
      source_cell <- strsplit(x = ccc.scores$ccc[ii], split = "=", fixed = TRUE)[[1]][1]
      target_cell <- strsplit(x = ccc.scores$ccc[ii], split = "=", fixed = TRUE)[[1]][2]
      
      lr_source <- temp2[[source_cell]]
      lr_target <- temp2[[target_cell]]
      
      ind1 <- which(lr_source$pfam_source=="PSEUDODOMAINLR")
      ind2 <- which(lr_target$pfam_source=="PSEUDODOMAINLR")
      
      if((length(ind1) > 0) && (length(ind2) > 0)){
        
        lr_source <- lr_source[ind1, ]
        lr_target <- lr_target[ind2, ]
        
        common_ligands <- intersect(x = c(lr_source$gene_source, lr_source$gene_target), 
                                    y = c(lr_target$gene_source, lr_target$gene_target))
        
        if(length(common_ligands) > 0){
          
          df_source <- temp2[[source_cell]]
          df_target <- temp2[[target_cell]]
          
          for(jj in 1:length(common_ligands)){
            
            idx <- which(df_source$gene_source==common_ligands[jj])
            if(length(idx) > 0){
              temp3[[source_cell]][idx] <- ccc.scores$score[ii]
            }
            idx <- which(df_source$gene_target==common_ligands[jj])
            if(length(idx) > 0){
              temp3[[source_cell]][idx] <- ccc.scores$score[ii]
            }
            
            idx <- which(df_target$gene_source==common_ligands[jj])
            if(length(idx) > 0){
              temp3[[target_cell]][idx] <- ccc.scores$score[ii]
            }
            idx <- which(df_target$gene_target==common_ligands[jj])
            if(length(idx) > 0){
              temp3[[target_cell]][idx] <- ccc.scores$score[ii]
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  calculate_weights <- function(lr, ccc) {
    if (is.na(lr) && is.na(ccc)) {
      return(NA)
    } else if (is.na(lr)) {
      return(ccc)
    } else if (is.na(ccc)) {
      return(lr)
    } else {
      return(lr * ccc)
    }
  }
  
  for(ii in 1:length(temp2)){
    temp2[[ii]]$ccc.scores <- temp3[[ii]]
    # temp2[[ii]]$lr.scores[which(is.na(temp2[[ii]]$lr.scores))] <- 1
    # temp2[[ii]]$ccc.scores[which(is.na(temp2[[ii]]$ccc.scores))] <- 1
    # temp2[[ii]]$weight <- temp2[[ii]]$lr.scores * temp2[[ii]]$ccc.scores
    temp2[[ii]]$weight <- mapply(calculate_weights, temp2[[ii]]$lr.scores, temp2[[ii]]$ccc.scores)
  }
  
  
  ### Now do the split between the joint receptors
  temp4 <- list()
  for(ii in 1:length(temp2)){
    
    curr <- temp2[[ii]]
    curr$lr.scores <- as.numeric(curr$lr.scores)
    curr$ccc.scores <- as.numeric(curr$ccc.scores)
    curr$weight <- as.numeric(curr$weight)
    receptors <- unique(curr$gene_target[which(grepl(pattern = "|", x = curr$gene_target, fixed = TRUE))])
    
    if(length(receptors)==0){
      
      temp4[[length(temp4)+1]] <- temp2[[ii]]
      
    } else {
      
      for(jj in 1:length(receptors)){
        
        proteins <- unique(unlist(strsplit(x = receptors[jj], split = "|", fixed = TRUE)))
        
        for(kk in 1:length(proteins)){
          
          domains <- unique(c(curr$pfam_source[which(curr$gene_source==proteins[kk])],
                              curr$pfam_target[which(curr$gene_target==proteins[kk])]))
          
          toBind <- matrix(data = , nrow = length(domains), ncol = 7)
          toBind[, 1] <- domains
          toBind[, 2] <- domains
          toBind[, 3] <- receptors[jj]
          toBind[, 4] <- proteins[kk]
          toBind[, 5] <- NA
          toBind[, 6] <- NA
          toBind[, 7] <- NA
          colnames(toBind) <- colnames(curr)
          toBind <- as.data.frame(toBind)
          toBind$lr.scores <- as.numeric(toBind$lr.scores)
          toBind$ccc.scores <- as.numeric(toBind$ccc.scores)
          toBind$weight <- as.numeric(toBind$weight)
          
          curr <- unique(rbind(curr, toBind))
            
        }
        
      }
      
      curr$lr.scores <- as.numeric(curr$lr.scores)
      curr$ccc.scores <- as.numeric(curr$ccc.scores)
      curr$weight <- as.numeric(curr$weight)
      temp4[[length(temp4)+1]] <- curr
      
    }
    
  }
  names(temp4) <- names(temp2)
  
  for(ii in 1:length(temp4)){
    temp4[[ii]]$lr.scores <- as.numeric(temp4[[ii]]$lr.scores)
    temp4[[ii]]$ccc.scores <- as.numeric(temp4[[ii]]$ccc.scores)
    temp4[[ii]]$weight <- as.numeric(temp4[[ii]]$weight)
  }
  background.networks.list$background.networks <- temp4
  
  return(background.networks.list)
  
}