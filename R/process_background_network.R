process_background_network <- function(background.networks.list = background.networks.list){
  
  bnList <- background.networks.list$background.networks
  temp <- list()
  for(ii in 1:length(bnList)){
    
    bn <- bnList[[ii]]
    idx <- intersect(x = which(is.na(bn$pfam_source)), 
                     y = which(is.na(bn$pfam_target)))
    
    ddi <- matrix(data = , nrow = 1, ncol = 4)
    for(jj in 1:length(idx)){
      
      domains <- c(bn$pfam_source[which(bn$gene_source==bn$gene_source[idx[jj]])],
                   bn$pfam_target[which(bn$gene_target==bn$gene_source[idx[jj]])])
      domains <- domains[complete.cases(domains)]
      
      if(length(domains) > 0){
        
        tobind <- matrix(data = , nrow = length(domains), ncol = 4)
        tobind[, 1] <- domains
        tobind[, 2] <- "DOMAIN"
        tobind[, 3] <- bn$gene_source[idx[jj]]
        tobind[, 4] <- bn$gene_target[idx[jj]]
        
        ddi <- rbind(ddi, tobind)
        
      }
      
    }
    
    if(nrow(ddi) > 1){
      ddi <- ddi[2:nrow(ddi), ]
      colnames(ddi) <- colnames(bn)
      ddi <- as.data.frame(ddi)
      bn <- rbind(bn, ddi)
    }
    
    idx2rem <- intersect(x = which(is.na(bn$pfam_source)), 
                         y = which(is.na(bn$pfam_target)))
    if(length(idx2rem) > 0){
      bn <- bn[-idx2rem, ]
    }
    
    idx <- which(is.na(bn$pfam_source))
    if(length(idx) > 0){
      bn$pfam_source[idx] <- "DOMAIN"
    }
    
    idx <- which(is.na(bn$pfam_target))
    if(length(idx) > 0){
      bn$pfam_target[idx] <- "DOMAIN"
    }
    
    bn <- bn[complete.cases(bn), ]
    
    temp[[length(temp)+1]] <- bn
    
  }
  names(temp) <- names(bnList)
  
  background.networks.list$background.networks <- temp
  
  return(background.networks.list)
  
}