#' Prepares the LINDAPlus object that contains information about the intra- and
#' inter-cellular interactions to train and specifies the list of receptors and
#' ligands.
#'
#'@param interactions_df 
#'
#'@param cell_types 
#'
#'@param genes2keep
#'
#'@param genes2remove
#'
#'@param receptors
#'
#'@param ligands
#'
#'
prepare_linda_inputs <- function(interactions_df = interactions_df,
                                 cell_types = cell_types,
                                 genes2keep = NULL,
                                 genes2remove = NULL,
                                 ligands = NULL,
                                 receptors = NULL){
  
  background.networks.list <- list()
  
  ## Ligand-Receptors
  ligands.receptors <- list()
  
  if((class(ligands) != "character") && (!is.null(ligands))){
    
    stop("Object 'cell_type' should be a character vector containing the names
         of each cell-type present in the multi-cellular system.")
    
  } else {
    
    if(is.null(ligands)){
      
      ind1 <- intersect(x = which(is.na(interactions_df$pfam_source)), 
                        y = which(is.na(interactions_df$pfam_target)))
      
      ind2 <- intersect(x = which(is.na(interactions_df$pfam_source)), 
                        y = which(!is.na(interactions_df$pfam_target)))
      
      ligands <- intersect(x = interactions_df$gene_target[ind1], 
                           y = interactions_df$gene_source[ind2])
      
      idx2rem1 <- which(!(interactions_df$gene_target[ind1]%in%ligands))
      idx2rem2 <- which(!(interactions_df$gene_source[ind2]%in%ligands))
      
      if(length(idx2rem1) > 0){
        interactions_df <- interactions_df[-ind1[idx2rem1]]
      }
      
      if(length(idx2rem2) > 0){
        interactions_df <- interactions_df[-ind2[idx2rem2]]
      }
      
    } else {
      
      ligands <- ligands
      
    }
    
  }
  
  if((class(receptors) != "character") && (!is.null(receptors))){
    
    stop("Object 'cell_type' should be a character vector containing the names
         of each cell-type present in the multi-cellular system.")
    
  } else {
    
    if(is.null(receptors)){
      
      ind <- intersect(x = which(is.na(interactions_df$pfam_source)), 
                       y = which(!is.na(interactions_df$pfam_target)))
      
      receptors <- unique(interactions_df$gene_target[ind])
      
    } else {
      
      receptors <- receptors
      
    }
    
  }
  
  idxdup2rem <- intersect(x = which(interactions_df$gene_source%in%receptors), 
                          y = which(interactions_df$gene_target%in%receptors))
  
  if(length(idxdup2rem) > 0){
    interactions_df <- interactions_df[-idxdup2rem, ]
  }
  
  ## Now prune ligands that are no TF targets
  ind1 <- intersect(x = which(is.na(interactions_df$pfam_source)), 
                    y = which(is.na(interactions_df$pfam_target)))
  ind2 <- intersect(x = which(is.na(interactions_df$pfam_source)), 
                    y = which(!is.na(interactions_df$pfam_target)))
  presLig <- unique(interactions_df$gene_target[ind1])
  lig2rem <- setdiff(x = presLig, y = ligands)
  if(length(lig2rem) > 0){
    interactions_df <- interactions_df[-which(interactions_df$gene_target%in%lig2rem), ]
  }
  presLig <- unique(interactions_df$gene_source[ind2])
  lig2rem <- setdiff(x = presLig, y = ligands)
  if(length(lig2rem) > 0){
    interactions_df <- interactions_df[-which(interactions_df$gene_source%in%lig2rem), ]
  }
  
  
  
  
  ## Prepare interactions for each cell-type
  if((class(cell_types) != "character") && (!is.null(cell_types))){
    stop("Object 'cell_type' should be a character vector containing the names
         of each cell-type present in the multi-cellular system.")
  }
  
  background.networks <- list()
  for(ii in 1:length(cell_types)){
    background.networks[[length(background.networks)+1]] <- interactions_df
  }
  names(background.networks) <- cell_types
  
  if(!is.null(genes2keep)){
    
    if((class(genes2keep) != "list") || 
       (is.null(names(genes2keep))) || 
       (length(genes2keep) != length(cell_types))){
      stop("Object 'genes2keep' should be either NULL (as set to default) or a 
           named list object for each cell-type containing the names of genes to 
           keep in the prior knowledge of interactions for each cell type.")
    } else {
      
      tmp_list <- list()
      nn <- names(genes2keep)
      for(ii in 1:length(nn)){
        
        curr <- background.networks[[nn[ii]]]
        idx2keep <- c(intersect(x = which(curr$gene_source%in%genes2keep[[ii]]), 
                                y = which(curr$gene_target%in%genes2keep[[ii]])),
                      which(is.na(curr$pfam_source)))
        tmp_idx_1 <- c()
        rec_idx <- which(grepl(pattern = "|", x = curr$gene_target, fixed = TRUE))
        if(length(rec_idx) > 0){
          for(jj in 1:length(rec_idx)){
            if(all(strsplit(x = curr$gene_target[rec_idx[jj]], split = "|", fixed = TRUE)[[1]]%in%genes2keep[[ii]])){
              tmp_idx_1 <- c(tmp_idx_1, rec_idx[jj])
            }
          }
        }
        tmp_idx_2 <- c()
        lig_idx <- which(grepl(pattern = "|", x = curr$gene_source, fixed = TRUE))
        if(length(lig_idx) > 0){
          for(jj in 1:length(lig_idx)){
            if(all(strsplit(x = curr$gene_source[lig_idx[jj]], split = "|", fixed = TRUE)[[1]]%in%genes2keep[[ii]])){
              tmp_idx_2 <- c(tmp_idx_2, lig_idx[jj])
            }
          }
        }
        idx2keep <- unique(c(idx2keep, tmp_idx_1, tmp_idx_2))
        if(length(idx2keep) > 0){
          
          curr <- curr[idx2keep, ]
          
        }
        tmp_list[[length(tmp_list)+1]] <- curr
        
      }
      
      background.networks <- tmp_list
      names(background.networks) <- nn
      
    }
    
  }
  
  if(!is.null(genes2remove)){
    
    if((class(genes2remove) != "list") || 
       (is.null(names(genes2remove))) || 
       (length(genes2remove) != length(cell_types))){
      stop("Object 'genes2remove' should be either NULL (as set to default) or a 
           named list object for each cell-type containing the names of genes to 
           keep in the prior knowledge of interactions for each cell type.")
    } else {
      
      tmp_list <- list()
      nn <- names(genes2remove)
      for(ii in 1:length(nn)){
        
        curr <- background.networks[[nn[ii]]]
        idx2rem <- c(which(curr$gene_source%in%genes2remove[[ii]]),
                     which(curr$gene_target%in%genes2remove[[ii]]))
        if(length(idx2rem) > 0){
          
          curr <- curr[-idx2rem, ]
          
        }
        tmp_list[[length(tmp_list)+1]] <- curr
        
      }
      
      background.networks <- tmp_list
      names(background.networks) <- nn
      
    }
    
  }
  
  tmp_list <- list()
  lig2keep <- c()
  for(ii in 1:length(background.networks)){
    
    bn <- background.networks[[ii]]
    
    # Process ligands
    ind1 <- intersect(x = which(is.na(bn$pfam_source)), 
                      y = which(!is.na(bn$pfam_target)))
    ind2 <- intersect(x = which(is.na(bn$pfam_source)), 
                      y = which(is.na(bn$pfam_target)))
    
    ll1 <- unique(bn$gene_source[ind1])
    ll2 <- unique(bn$gene_target[ind2])
    
    ligands <- intersect(x = ll1, y = ll2)
    
    lig2rem <- setdiff(x = ll1, y = ligands)
    if(length(lig2rem) > 0){
      bn <- bn[-ind1[which(bn$gene_source[ind1]%in%lig2rem)], ]
    }
    
    lig2rem <- setdiff(x = ll2, y = ligands)
    if(length(lig2rem) > 0){
      bn <- bn[-ind2[which(bn$gene_target[ind2]%in%lig2rem)], ]
    }
    
    lig2keep <- c(lig2keep, ligands)
    
    tmp_list[[length(tmp_list)+1]] <- bn
    
  }
  names(tmp_list) <- names(background.networks)
  ligands <- unique(lig2keep)
  
  rec2rem <- c()
  split_receptors <- unique(unlist(strsplit(x = receptors, split = "|", fixed = TRUE)))
  for(ii in 1:length(background.networks)){
    
    curr <- background.networks[[ii]]
    curr_rec <- unique(curr$gene_target[intersect(x = which(is.na(curr$pfam_source)), 
                                                  y = which(!is.na(curr$pfam_target)))])
    curr_rec_split <- unique(unlist(strsplit(x = curr_rec, split = "|", fixed = TRUE)))
    rec2rem <- c(rec2rem, setdiff(x = split_receptors, y = curr_rec_split))
    
  }
  
  bnList <- list()
  for(ii in 1:length(background.networks)){
    
    curr <- background.networks[[ii]]
    idx2rem <- which(curr$gene_source%in%rec2rem)
    if(length(idx2rem) > 0){curr <- curr[-idx2rem, ]}
    bnList[[length(bnList)+1]] <- curr
    
  }
  names(bnList) <- names(background.networks)
  background.networks <- bnList
  
  rec2keep <- setdiff(x = split_receptors, y = rec2rem)
  rec2keep1 <- receptors[which(receptors%in%rec2keep)]
  rec2keep2 <- c()
  for(ii in 1:length(rec2keep)){
    rec2keep2 <- c(rec2keep2, receptors[which(grepl(pattern = paste0(rec2keep[ii], "|"), x = receptors, fixed = TRUE))])
  }
  rec2keep3 <- c()
  for(ii in 1:length(rec2keep)){
    rec2keep3 <- c(rec2keep3, receptors[which(grepl(pattern = paste0("|", rec2keep[ii]), x = receptors, fixed = TRUE))])
  }
  
  receptors <- unique(c(rec2keep1, rec2keep2, rec2keep3))
  
  ligands.receptors[[length(ligands.receptors)+1]] <- ligands
  ligands.receptors[[length(ligands.receptors)+1]] <- receptors
  names(ligands.receptors) <- c("ligands", "receptors")
  
  background.networks.list[[length(background.networks.list)+1]] <- background.networks
  background.networks.list[[length(background.networks.list)+1]] <- ligands.receptors
  names(background.networks.list) <- c("background.networks", "ligands.receptors")
  
  return(background.networks.list)
  
}
