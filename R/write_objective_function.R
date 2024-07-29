write_objective_function <- function(variables = variables,
                                     background.networks.list = background.networks.list,
                                     tf.scores = tf.scores,
                                     ligand.scores = ligand.scores,
                                     lr.scores = lr.scores,
                                     ccc.scores = ccc.scores,
                                     lambda1 = lambda1, 
                                     lambda2 = lambda2, 
                                     lambda3 = lambda3,
                                     lambda4 = lambda4){
  
  print("Writing the objective function and constraints. This might take a bit of time..")
  
  cells <- names(tf.scores)
  background.network.list <- background.networks.list$background.networks
  
  of1 <- ""
  for(kk in 1:length(cells)){
    
    # Write main objective - Minimizing the TF scores
    input.scores <- tf.scores[[kk]]
    tfNodes <- input.scores$tf
    idx <- which(variables$var_exp%in%paste0(cells[kk], ":node ", tfNodes))
    if(length(idx)==0){
      stop(paste0("No input TF mapped to the background network for cell-type ", cells[kk], ". ",
         "Check inputs/background network"))
    }
    
    objective.function <- ""
    cnt <- 1
    for(ii in seq_len(nrow(input.scores))){
      
      idx <- which(variables$var_exp==paste0(cells[kk], ":node ", input.scores$tf[ii]))
      if(length(idx)==1){
        
        if(cnt == 1){
          
          if(input.scores$score[ii]==0){
            objective.function <- paste0(objective.function,
                                         lambda1,
                                         " ",
                                         variables$var[idx])
          } else {
            objective.function <- paste0(objective.function,
                                         "- ",
                                         lambda1,
                                         " ",
                                         variables$var[idx])
          }
          cnt <- cnt + 1
          
        } else {
          
          if(input.scores$score[ii]==0){
            objective.function <- paste0(objective.function,
                                         " + ",
                                         lambda1,
                                         " ",
                                         variables$var[idx])
          } else {
            objective.function <- paste0(objective.function,
                                         " - ",
                                         lambda1,
                                         " ",
                                         variables$var[idx])
          }
          
        }
        
      }
      
    }
    
    of1 <- paste0(of1, objective.function, collapse = " ")
    
  }
  
  # The secondary Ligands Objective Function
  if(is.null(ligand.scores)){
    
    of2 <- ""
    
  } else {
    
    topLigands <- ligand.scores$ligand[which(ligand.scores$score==1)]
    penLigands <- ligand.scores$ligand[which(ligand.scores$score==0)]
    of2 <- ""
    ofLigands <- c(topLigands, penLigands)
    for(ii in 1:length(ofLigands)){
      
      if(ofLigands[ii]%in%topLigands){
        of2 <- paste0(of2, " - ", lambda2, " ", variables$var[which(variables$var_exp==paste0("LR:ligand ", ofLigands[ii]))])
      }
      
      if(ofLigands[ii]%in%penLigands){
        of2 <- paste0(of2, " + ", lambda2, " ", variables$var[which(variables$var_exp==paste0("LR:ligand ", ofLigands[ii]))])
      }
      
    }
    
  }
  objective.function <- paste0(of1, of2, collapse = " + ")
  
  if((is.null(lr.scores)) && (is.null(ccc.scores))){
    
    of3 <- ""
    
  } else {
    
    of3 <- ""
    cc <- c()
    for(ii in 1:length(background.network.list)){
      
      df <- background.network.list[[ii]]
      
      ind <- intersect(x = intersect(x = which(df$pfam_source == "PSEUDODOMAINLR"), 
                                     y = which(df$pfam_target != "PSEUDODOMAINLR")), 
                       y = which(!is.na(df$weight)))
      
      if(length(ind) > 0){
        
        ppp <- unique(paste0(df$gene_source[ind], "=", df$gene_target[ind]))
        ppp_src <- paste0(df$gene_source, "=", df$gene_target)
        
        varvar1 <- rep("", length(ppp))
        valval1 <- rep(0, length(ppp))
        for(jj in 1:length(ppp)){
          
          varvar1[jj] <- variables$var[which(variables$var_exp == paste0(names(background.network.list)[ii], 
                                                                         ":interaction ", 
                                                                         strsplit(ppp[jj], split = "=", fixed = TRUE)[[1]][1], 
                                                                         "=", 
                                                                         strsplit(ppp[jj], split = "=", fixed = TRUE)[[1]][2]))]
          
          if(all(is.na(df$weight[which(ppp_src==ppp[jj])]))){
            valval1[jj] <- 0
          } else {
            valval1[jj] <- mean(df$weight[which(ppp_src==ppp[jj])], na.rm = TRUE)
          }
          
        }
        
        fff <- rep("", length(ppp))
        for(jj in 1:length(ppp)){
          fff[jj] <- paste0(lambda3 * abs(valval1[jj]), " ", varvar1[jj])
        }
        of3 <- paste0(fff, collapse = " - ")
        of3 <- paste0(" - ", of3)
        
        cc <- c(cc, of3)
        
        
      }
      
    }
    
    objective.function <- paste0(objective.function, paste0(cc, collapse = ""))
    
  }
  
  
  
  # Write fourth objective - size penalty factor
  
  # idx <- c(which(grepl(pattern = "interaction", x = variables$var_exp, fixed = TRUE)),
  #          intersect(x = which(grepl(pattern = "reaction ", x = variables$var_exp, fixed = TRUE)),
  #                    y = which(grepl(pattern = "PSEUDODOMAIN", x = variables$var_exp, fixed = TRUE))))
  # idx <- which(grepl(pattern = "interaction", x = variables$var_exp, fixed = TRUE))
  
  # idx <- which(grepl(pattern = "reaction", x = variables$var_exp, fixed = TRUE))
  
  idx <- which(grepl(pattern = "domain", x = variables$var_exp, fixed = TRUE))
  
  if(lambda4>0){
    obj <- paste0(" + ", lambda4, " ", variables$var[idx])
    obj <- paste0(obj, collapse = "")
    objective.function <- paste0(objective.function, obj)
  } else {
    if(lambda4<0){
      obj <- paste0(" - ", abs(lambda4), " ", variables$var[idx])
      obj <- paste0(obj, collapse = "")
      objective.function <- paste0(objective.function, obj)
    }
  }
  
  return(objective.function)
  
}