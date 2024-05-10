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
        idx2keep <- intersect(x = which(curr$gene_source%in%genes2keep[[ii]]), 
                              y = which(curr$gene_target%in%genes2keep[[ii]]))
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
  
  
  ## Ligand-Receptors
  ligands.receptors <- list()
  
  if((class(ligands) != "character") && (!is.null(ligands))){
    
    stop("Object 'cell_type' should be a character vector containing the names
         of each cell-type present in the multi-cellular system.")
    
  } else {
    
    if(is.null(ligands)){
      
      ind <- intersect(x = which(is.na(interactions_df$pfam_source)), 
                       y = which(!is.na(interactions_df$pfam_target)))
      
      ligands <- unique(interactions_df$gene_source[ind])
      
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
      
    }
    
  }
  
  ligands.receptors[[length(ligands.receptors)+1]] <- ligands
  ligands.receptors[[length(ligands.receptors)+1]] <- receptors
  names(ligands.receptors) <- c("ligands", "receptors")
  
  background.networks.list[[length(background.networks.list)+1]] <- background.networks
  background.networks.list[[length(background.networks.list)+1]] <- ligands.receptors
  names(background.networks.list) <- c("background.networks", "ligands.receptors")
  
  return(background.networks.list)
  
}