write_objective_function <- function(variables = variables, 
                                     background.networks.list = background.networks.list, 
                                     ccc.input = ccc.input,
                                     tf.input = tf.input,
                                     ligand.scores = ligand.scores,
                                     ccc.prob = ccc.prob,
                                     lambda1 = lambda1, 
                                     lambda2 = lambda2, 
                                     lambda3 = lambda3,
                                     lambda4 = lambda4,
                                     lambda5 = lambda5){
  
  print("Writing the objective function and constraints. This might take a bit of time..")
  
  cells <- names(background.networks.list$background.networks)
  background.network.list <- background.networks.list$background.networks
  
  
  if(is.null(ccc.input)){
    
    of1 <- ""
    of2 <- ""
    
  } else {
    
    #### Write main objective function: Minimize ligands in transmitter cells
    of1 <- ""
    transmitter_ligands <- c()
    for(ii in 1:length(cells)){
      transmitter_ligands <- c(transmitter_ligands, intersect(x = variables$var_exp, 
                                                              y = paste0(cells[ii], ":node ", background.networks.list$ligands.receptors$ligands)))
    }
    
    transmitter_vars <- rep("", length(transmitter_ligands))
    for(ii in 1:length(transmitter_ligands)){
      transmitter_vars[ii] <- variables$var[which(variables$var_exp == transmitter_ligands[ii])]
    }
    
    transmitter_vals <- rep(lambda1, length(transmitter_vars))
    
    for(ii in 1:nrow(ccc.input)){
      ind <- which(transmitter_ligands == paste0(ccc.input$transmitter_cell[ii], ":node ", ccc.input$ligand[ii]))
      if(length(ind) > 0){
        transmitter_vals[ind] <- (-1)*lambda1
      }
    }
    
    for(ii in 1:length(transmitter_vals)){
      
      if(ii == 1){
        
        of1 <- paste0(transmitter_vals[ii], " ", transmitter_vars[ii])
        
      } else {
        
        if(transmitter_vals[ii] < 0){
          
          of1 <- paste0(of1, " - ", abs(transmitter_vals[ii]), " ", transmitter_vars[ii])
          
        } else {
          
          of1 <- paste0(of1, " + ", abs(transmitter_vals[ii]), " ", transmitter_vars[ii])
          
        }
        
      }
      
    }
    
    
    
    #### Secondary objective function: ligand receptor interaction based on ccc.input evidence
    of2 <- ""
    lr_int <- variables$var_exp[which(grepl(pattern = "LR:interaction ", x = variables$var_exp, fixed = TRUE))]
    lr_int <- unique(gsub(pattern = "LR:interaction ", replacement = "", x = lr_int, fixed = TRUE))
    lr_exp <- c()
    lr_var <- c()
    for(ii in 1:length(lr_int)){
      for(jj in 1:length(cells)){
        ind <- which(variables$var_exp == paste0(cells[jj], ":interaction ", lr_int[ii]))
        if(length(ind) > 0){
          lr_exp <- c(lr_exp, variables$var_exp[ind])
          lr_var <- c(lr_var, variables$var[ind])
        }
      }
    }
    lr_val <- rep(lambda2, length(lr_var))
    
    for(ii in 1:nrow(ccc.input)){
      ind <- which(lr_exp == paste0(ccc.input$receiver_cell[ii], ":interaction ", ccc.input$ligand[ii], "=", ccc.input$receptor[ii]))
      if(length(ind) > 0){
        lr_val[ind] <- (-1)*lambda2
      }
    }
    
    for(ii in 1:length(lr_val)){
      
      if(ii == 1){
        
        of2 <- paste0(lr_val[ii], " ", lr_var[ii])
        
      } else {
        
        if(lr_val[ii] < 0){
          
          of2 <- paste0(of2, " - ", abs(lr_val[ii]), " ", lr_var[ii])
          
        } else {
          
          of2 <- paste0(of2, " + ", abs(lr_val[ii]), " ", lr_var[ii])
          
        }
        
      }
      
    }
    
  }
  
  
  
  #### Third objective function: TF scores
  if(is.null(tf.input)){
    
    of4 <- ""
    
  } else {
    
    cells <- names(tf.input)
    background.network.list <- background.networks.list$background.networks
    
    of4 <- ""
    for(kk in 1:length(cells)){
      
      # Write main objective - Minimizing the TF scores
      input.scores <- tf.input[[kk]]
      tfNodes <- input.scores$tf
      idx <- which(variables$var_exp%in%paste0(cells[kk], ":node ", tfNodes))
      if(length(idx)==0){
        stop(paste0("No input TF mapped to the background network for cell-type ", cells[kk], ". ",
                    "Check inputs/background network"))
      }
      
      of <- ""
      cnt <- 1
      for(ii in seq_len(nrow(input.scores))){
        
        idx <- which(variables$var_exp==paste0(cells[kk], ":node ", input.scores$tf[ii]))
        if(length(idx)==1){
          
          if(cnt == 1){
            
            if(input.scores$score[ii]==0){
              of <- paste0(of, lambda3, " ", variables$var[idx])
            } else {
              of <- paste0(of, "- ", lambda3, " ", variables$var[idx])
            }
            cnt <- cnt + 1
            
          } else {
            
            if(input.scores$score[ii]==0){
              of <- paste0(of, " + ", lambda3, " ", variables$var[idx])
            } else {
              of <- paste0(of, " - ", lambda3, " ", variables$var[idx])
            }
            
          }
          
        }
        
      }
      
      of4 <- paste0(of4, " ", of)
      
    }
    
    
  }
  
  
  
  
  #### Third objective function: Ligand scores
  if(is.null(ligand.scores)){
    
    of3 <- ""
    
  } else {
    
    of3 <- ""
    cnt <- 1
    for(ii in 1:nrow(ligand.scores)){
      
      ind <- which(variables$var_exp == paste0("LR:ligand ", ligand.scores$ligand[ii]))
      if(length(ind) > 0){
        if(cnt == 1){
          of3 <- paste0(of3, lambda4*(1-ligand.scores$score[ii]), " ", variables$var[ind])
          cnt <- cnt + 1
        } else {
          of3 <- paste0(of3, " + ", lambda4*(1-ligand.scores$score[ii]), " ", variables$var[ind])
        }
      }
    }
  }
  
  
  
  
  
  
  #### Write last objective - size penalty factor
  if((of1 == "") && (of4 != "")){
    
    objective.function <- substr(x = of4, start = 2, stop = nchar(of4))
    
  } else {
    
    if((of1 != "") && (of4 == "")){
      
      objective.function <- of1
      if(substr(x = of2, start = 1, stop = 1) == "-"){
        objective.function <- paste0(objective.function, " - ", substr(x = of2, start = 2, stop = nchar(of2)))
      } else {
        objective.function <- paste0(objective.function, " + ", of2)
      }
      
    } else {
      
      objective.function <- of1
      if(substr(x = of2, start = 1, stop = 1) == "-"){
        objective.function <- paste0(objective.function, " - ", substr(x = of2, start = 2, stop = nchar(of2)))
      } else {
        objective.function <- paste0(objective.function, " + ", of2)
      }
      
      if(substr(x = of4, start = 2, stop = 2) == "-"){
        objective.function <- paste0(objective.function, " ", of4)
      } else {
        objective.function <- paste0(objective.function, " + ", of4)
      }
      
    }
    
  }
  objective.function <- gsub(pattern = "  ", replacement = " ", fixed = TRUE, x = objective.function)
  
  if(of3 != ""){
    objective.function <- paste0(objective.function, " + ", of3)
  }
  
  
  idx <- which(grepl(pattern = "domain", x = variables$var_exp, fixed = TRUE))
  
  if(lambda5>0){
    obj <- paste0(" + ", lambda5, " ", variables$var[idx])
    obj <- paste0(obj, collapse = "")
    objective.function <- paste0(objective.function, obj)
  } else {
    if(lambda5<0){
      obj <- paste0(" - ", abs(lambda5), " ", variables$var[idx])
      obj <- paste0(obj, collapse = "")
      objective.function <- paste0(objective.function, obj)
    }
  }
  
  return(objective.function)
  
}