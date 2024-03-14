write_objective_function <- function(variables = variables,
                                     background.networks.list = background.networks.list,
                                     tf_scores = tf_scores,
                                     ligand_scores = ligand_scores,
                                     alpha = 10,
                                     beta = 5,
                                     gamma = 0.1){
  
  print("Writing the objective function and constraints.
        This might take a bit of time..")
  
  cells <- names(tf_scores)
  
  of1 <- ""
  for(kk in 1:length(cells)){
    
    # Write main objective - Minimizing the TF scores
    input.scores <- tf_scores[[kk]]
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
                                         alpha,
                                         " ",
                                         variables$var[idx])
          } else {
            objective.function <- paste0(objective.function,
                                         "- ",
                                         alpha,
                                         " ",
                                         variables$var[idx])
          }
          cnt <- cnt + 1
          
        } else {
          
          if(input.scores$score[ii]==0){
            objective.function <- paste0(objective.function,
                                         " + ",
                                         alpha,
                                         " ",
                                         variables$var[idx])
          } else {
            objective.function <- paste0(objective.function,
                                         " - ",
                                         alpha,
                                         " ",
                                         variables$var[idx])
          }
          
        }
        
      }
      
    }
    
    of1 <- paste0(of1, objective.function, collapse = " ")
    
  }
  
  # The secondary Ligands Objective Function
  topLigands <- ligand_scores$ligand[which(ligand_scores$score==1)]
  penLigands <- ligand_scores$ligand[which(ligand_scores$score==0)]
  of2 <- ""
  ofLigands <- c(topLigands, penLigands)
  for(ii in 1:length(ofLigands)){
    
    if(ofLigands[ii]%in%topLigands){
      of2 <- paste0(of2, " - ", beta, " ", variables$var[which(variables$var_exp==paste0("LR:ligand ", ofLigands[ii]))])
    }
    
    if(ofLigands[ii]%in%penLigands){
      of2 <- paste0(of2, " + ", beta, " ", variables$var[which(variables$var_exp==paste0("LR:ligand ", ofLigands[ii]))])
    }
    
  }
  objective.function <- paste0(of1, of2)
  
  # Write third objective - size penalty factor
  idx <- which(grepl(pattern = ":domain ", x = variables$var_exp, fixed = TRUE))
  if(gamma>0){
    obj <- paste0(" + ", gamma, " ", variables$var[idx])
    obj <- paste0(obj, collapse = "")
    objective.function <- paste0(objective.function, obj)
  } else {
    if(gamma<0){
      obj <- paste0(" - ", abs(gamma), " ", variables$var[idx])
      obj <- paste0(obj, collapse = "")
      objective.function <- paste0(objective.function, obj)
    }
  }
  
  return(objective.function)
  
}