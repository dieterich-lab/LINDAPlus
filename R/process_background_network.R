process_background_network <- function(background.networks.list = background.networks.list, 
                                       tf.input = tf.input,
                                       ccc.prob = ccc.prob,
                                       ccc.input = ccc.input){
  
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
    
    idx <- intersect(x = which(bn$gene_source%in%tf.input[[ii]]$tf), 
                     y = which(bn$gene_target%in%background.networks.list$ligands.receptors$ligands))
    
    if(length(idx) > 0){
      bn$pfam_source[idx] <- "PSEUDODOMAINTF"
      bn$pfam_target[idx] <- "PSEUDODOMAINTF"
    }
    
    temp[[length(temp)+1]] <- bn
    
  }
  names(temp) <- names(bnList)
  
  temp2 <- list()
  for(ii in 1:length(temp)){
    
    curr <- temp[[ii]]
    curr$lr.scores <- NA
    curr$ccc.scores <- NA
    curr$weight <- NA
    temp2[[length(temp2)+1]] <- curr
    
  }
  names(temp2) <- names(bnList)
  
  temp3 <- list()
  for(ii in 1:length(temp2)){
    
    curr <- temp2[[ii]]
    ind <- which(curr$pfam_target == "PSEUDODOMAINTF")
    if(length(ind) > 0){
      
      curr$pfam_source[ind] <- "PSEUDODOMAINTF"
      
    }
    temp3[[length(temp3)+1]] <- unique(curr)
    
  }
  names(temp3) <- names(bnList)
  
  background.networks.list$background.networks <- temp3
  
  return(background.networks.list)
  
}