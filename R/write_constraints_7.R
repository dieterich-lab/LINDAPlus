write_constraints_7 <- function(variables = variables, 
                                background.networks.list = background.networks.list,
                                as.input = as.input){
  
  if(is.null(as.input)){
    
    constraints <- c()

  } else {
    
    ## Exclusion
    cc1 <- c()
    idx <- which(as.input$effect=="exclusion")
    if(length(idx) > 0){
      
      for(ii in 1:length(idx)){
        
        # domain
        ind <- which(variables$var_exp==paste0(as.input$cell_type[idx[ii]],
                                               ":domain ",
                                               as.input$domainID[idx[ii]],
                                               " of protein ",
                                               as.input$proteinID[idx[ii]]))
        cc1 <- c(cc1, paste0(variables$var[ind], " = 0"))
        
      }
      
    }
    
    
    ## Inclusion
    cc2 <- c()
    idx <- which(as.input$effect=="inclusion")
    if(length(idx) > 0){
      
      for(ii in 1:length(idx)){
        
        bn <- background.networks.list$background.networks[[as.input$cell_type[idx[ii]]]]
        indss <- intersect(x = which(bn$pfam_source==as.input$domainID[idx[ii]]), 
                           y = which(bn$gene_source==as.input$proteinID[idx[ii]]))
        indtt <- intersect(x = which(bn$pfam_target==as.input$domainID[idx[ii]]), 
                           y = which(bn$gene_target==as.input$proteinID[idx[ii]]))
        
        if(length(indss) > 0){
          
          tmp <- bn[indss, ]
          
        }
        
        if(length(indtt) > 0){
          
          tmp <- bn[indtt, ]
          idx2rem <- c()
          pp <- paste0(as.input$cell_type, "::", as.input$proteinID, "::", 
                       as.input$domainID, "::", as.input$effect)
          for(jj in 1:nrow(tmp)){
            
            nn <- intersect(x = pp, 
                            y = paste0(as.input$cell_type[idx[ii]], "::",
                                       tmp$gene_source[jj], "::",
                                       tmp$pfam_source[jj], "::exclusion"))
            if(length(nn) > 0){
              idx2rem <- c(idx2rem, jj)
            }
            
          }
          
          if(length(idx2rem) < nrow(tmp)){
            
            vv <- c()
            if(length(idx2rem) > 0){
              
              tmp <- tmp[-idx2rem, ]
              
            }
            
            for(jj in 1:nrow(tmp)){
              
              ind <- which(variables$var_exp==paste0(as.input$cell_type[idx[ii]],
                                                     ":reaction ",
                                                     tmp$pfam_source[jj],
                                                     "=",
                                                     tmp$pfam_target[jj],
                                                     " of ",
                                                     tmp$gene_source[jj],
                                                     "=",
                                                     tmp$gene_target))
              vv <- c(vv, variables$var[ind])
              
            }
            
            cc2 <- c(cc2, paste0(paste0(vv, collapse = " + "), " > 0"))
            
            cc2 <- c(cc2, paste0(variables$var[which(variables$var_exp==
                                                       paste0(as.input$cell_type[idx[ii]],
                                                              ":domain ",
                                                              as.input$domainID[idx[ii]],
                                                              " of protein ",
                                                              as.input$proteinID[idx[ii]]))],
                                 " = 1"))
            
          }
          
          
        }
        
        
      }
      
    }
    
    constraints <- unique(c(cc1, cc2))
    
  }
  
  return(constraints)
  
}