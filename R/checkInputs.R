checkInputs <- function(background.networks.list = background.networks.list,
                        ccc.input = ccc.input,
                        tf.input = tf.input,
                        top.tf = top.tf,
                        ligand.scores = ligand.scores,
                        ccc.prob = ccc.prob,
                        as.input = as.input,
                        solverPath = solverPath,
                        lambda1 = lambda1,
                        lambda2 = lambda2,
                        lambda3 = lambda3,
                        lambda4 = lambda4,
                        lambda5 = lambda5,
                        mipgap = mipgap,
                        relgap = relgap,
                        populate = populate,
                        nSolutions = nSolutions,
                        timelimit = timelimit,
                        intensity = intensity,
                        replace = replace,
                        threads = threads,
                        condition = condition,
                        save_res = save_res){ 
  
  #### background.networks.list
  if((class(background.networks.list) != "list") ||
     (length(background.networks.list)!= 2) || 
     (length(intersect(x = names(background.networks.list), y = c("background.networks", "ligands.receptors")))<2)){
    
    stop("Object 'background.networks.list' should be a list with a length of 2 named 'background.networks' and 'ligands.receptors'. Please check your inputs.")
    
  } else {
    
    obj1 <- background.networks.list$background.networks
    obj2 <- background.networks.list$ligands.receptors
    
    #obj1
    cc <- c()
    for(ii in 1:length(obj1)){cc <- c(cc, class(obj1[[ii]]))}
    if((is.null(names(obj1))) ||
       all(cc != "data.frame")){
      stop("The 'background.networks' should be named data-frames. Please check your inputs.")
    } else {
      
      for(ii in 1:length(obj1)){
        
        df <- obj1[[ii]]
        if(length(intersect(x = colnames(df), y = c("pfam_source", "pfam_target", "gene_source", "gene_target")))<4){
          
          stop("The 'background.network' object should be a data-frame with at least 4 columns and have at lease c('pfam_source', 'pfam_target', 'gene_source', 'gene_target') as column ID's. Please check your inputs.")
        }
        
      }
    }
    
    #obj2
    cc <- c()
    for(ii in 1:length(obj2)){cc <- c(cc, class(obj2[[ii]]))}
    if(length(intersect(x = names(obj2), y = c("ligands", "receptors")))<2 ||
       all(cc != "character")){
      stop("The 'ligands.receptors' should be a list of characters with a length 
           of at least 2 named 'ligands' and 'receptors' and optionally 
           'perturbation.ligands'. Please check your inputs.")
    }
    
  }
  
  
  if(is.null(ccc.input) && is.null(tf.input)){
    
    stop("You should either provide a table of lignand-receptor interaction of
         cell-cell commuincation or a list of tables of TF score enrichments
         for each cell-type. Please check your inputs or refer to the tutorials
         about how to format your inputs.")
    
  }
  
  
  ### ccc.input
  if(!is.null(ccc.input)){
    
    if(class(ccc.input) != "data.frame"){
      
      stop("The 'ccc.input' object should be a data-frame with column names:
         'transmitter_cell', 'receiver_cell', 'ligand', 'receptor'. Please check your inputs.")
      
    } else {
      
      cn <- intersect(x = colnames(ccc.input), y = c("transmitter_cell", "receiver_cell", "ligand", "receptor"))
      if(length(cn) < 4){
        
        stop("The 'ccc.input' object should be a data-frame with column names:
         'transmitter_cell', 'receiver_cell', 'ligand', 'receptor'. Please check your inputs.")
        
      } else {
        
        nn1 <- setdiff(x = ccc.input$transmitter_cell, y = names(background.networks.list$background.networks))
        nn2 <- setdiff(x = ccc.input$receiver_cell, y = names(background.networks.list$background.networks))
        
        if(length(nn1) > 0){
          
          stop(paste0("The ", paste0(nn1, collapse = ", "), " cell-types of the transmitter cell have not been defined in the 'background.networks.list' object. ",
                      "Please check your inputs."))
          
        } else {
          
          if(length(nn2) > 0){
            
            stop(paste0("The ", paste0(nn2, collapse = ", "), " cell-types of the receiver cell have not been defined in the 'background.networks.list' object. ",
                        "Please check your inputs."))
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  #### tf.input
  if(!is.null(tf.input)){
    
    if((class(tf.input) != "list") ||
       (length(tf.input)!= length(background.networks.list$background.networks)) || 
       (length(intersect(x = names(background.networks.list$background.networks), y = names(tf.input)))!=length(tf.input))){
      
      stop("The 'tf.input' object should be a list of data-frames with the same length and named the same as the 'background.networks.list$background.networks' object. Please check your inputs.")
      
    } else {
      
      cc <- c()
      for(ii in 1:length(tf.input)){cc <- c(cc, class(tf.input[[ii]]))}
      if(all(cc != "data.frame")){
        stop("The 'tf.input' object should be a list of data-frames with the same length and named the same as the 'background.networks.list$background.networks' object. Please check your inputs.")
      }
      
    }
    
    if(is.null(top.tf)){
      
      warning("You have provided the 'tf.input' object, but you did not provide 
              the number of significant TF's for each cell-type in 'top.tf'. All 
              provided TF's for each cell-type will be considered as 
              significant.")
      
      top.tf <- rep(1, length(tf.input))
      for(ii in 1:length(tf.input)){
        top.tf[ii] <- nrow(tf.input[[ii]])
      }
      
      names(top.tf) <- names(tf.input)
      
    } else {
      
      if((!is.numeric(top.tf)) ||
         (length(intersect(x = names(top.tf), y = names(background.networks.list$background.networks)))<length(top.tf))){
        
        stop("The 'top.tf' parameter should be a named numeric vector with 
             cell-types given as names. Please check your inputs.")
        
      } else {
        
        temp <- list()
        for(ii in 1:length(top.tf)){
          
          curr <- tf.input[[which(names(tf.input)==names(top.tf)[ii])]]
          
          if(top.tf[ii] > nrow(curr)){
            
            warning(paste0("There a re more given top TF's for cell-type ", 
                           names(top.tf)[ii], " than there are TF's given. All 
                           the given TF's will be considered as significant."))
            
            curr$score <- 1
            temp[[length(temp)+1]] <- curr
            
          } else {
            
            curr <- curr[order(abs(as.numeric(curr$score)), decreasing = TRUE), ]
            curr$score <- 0
            curr$score[1:top.tf[ii]] <- 1
            temp[[length(temp)+1]] <- curr
            
          }
          
        }
        names(temp) <- names(top.tf)
        tf.input <- temp
        
      }
      
    }
    
  }
  
  
  #### ligand.scores
  if(!is.null(ligand.scores)){
    
    if(class(ligand.scores) != "data.frame"){
      
      stop("The 'ligand.scores' should be a data-frame object with column names 'ligand' and 'score' where the score values
           should be normalized between 0 and 1. Please check your inputs.")
      
    } else {
      
      nn <- intersect(x = colnames(ligand.scores), y = c("ligand", "score"))
      if(length(nn) < 2){
        
        stop("The 'ligand.scores' should be a data-frame object with column names 'ligand' and 'score' where the score values
           should be normalized between 0 and 1. Please check your inputs.")
        
      } else {
        
        if(!all(ligand.scores$score >= 0 & ligand.scores$score <= 1)){
          
          stop("The 'ligand.scores' should be a data-frame object with column names 'ligand' and 'score' where the score values
           should be normalized between 0 and 1. Please check your inputs.")
          
        } else {
          
          nn <- setdiff(x = unique(ligand.scores$ligand), y = background.networks.list$ligands.receptors$ligands)
          if(length(nn) == nrow(ligand.scores)){
            
            stop("None of the ligands specified in the 'ligand.scores' object are present in the 'background.networks.list$ligands.receptors$ligands'.
                 Please check your inputs.")
            
          } else {
            
            if(length(nn) > 0){
              
              warning(paste0("Ligands ", paste0(nn, collapse = ", "), " of the 'ligand.scores' object are not present in the
                             'background.networks.list$ligands.receptors$ligands' and shall be removed."))
              
              ind <- which(ligand.scores$ligand %in% nn)
              ligand.scores <- ligand.scores[-ind, ]
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  #### ccc.prob
  if(!is.null(ccc.prob)){
    
    if(class(ccc.prob) != "data.frame"){
      
      stop("The provided 'ccc.prob' object should be either a data-frame with colnames 'ccc' (character) and 'scores' (numerical) or set to NULL (default). Please check on your inputs.")
      
    } else {
      
      cc <- intersect(x = colnames(ccc.prob), y = c("ccc", "score"))
      if(length(cc) < 2){
        
        stop("The provided 'ccc.prob' object should be either a data-frame with colnames 'ccc' (character) and 'scores' (numerical) or set to NULL (default). Please check on your inputs.")
        
      } else {
        
        if((class(ccc.prob$ccc) != "character") || 
           (class(ccc.prob$score) != "numeric")){
          
          stop("The provided 'ccc.prob' object should be either a data-frame with colnames 'ccc' (character) and 'scores' (numerical) or set to NULL (default). Please check on your inputs.")
          
        } else {
          
          cell_types <- unique(unlist(strsplit(x = ccc.prob$ccc, split = "=", fixed = TRUE)))
          scores <- ccc.prob$score
          
          cc <- intersect(x = cell_types, y = names(background.networks.list$background.networks))
          if(length(cc) < length(background.networks.list$background.networks)){
            
            stop("You should provide cell-cell communication probability scores for each pair of cell-types present in the system separated by a '=' symbol (i.e. 'CellA=CellB'). Please check your inputs.")
            
          }
          
          if(!((all(scores >= 0)) && all(scores <= 1))){
            
            stop("You should provide cell-cell communication probability scores for each pair of cell-types present in the system as numerical values between 0 and 1. Please check your inputs.")
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  #### as.input
  if(!is.null(as.input)){
    
    if((class(as.input) != "data.frame") ||
       (ncol(as.input) < 4) ||
       (length(intersect(x = colnames(as.input), 
                         y = c("cell_type", "proteinID", "domainID", "effect")))<4)){
      
      stop("If you use the 'as.input' input, it should be provided as a data-frame with 4 columns a colnames: 'cell_type', 'proteinID', 'domainID' and 'effect'. Otherwise, you can set 'as.input=NULL' in order to not account for splicing effects.")
      
    } else {
      
      if(!all(as.input$effect%in%c("exclusion", "inclusion"))){
        
        stop("In the 'effect' column of the 'as.input' data-frame object, users should either give character values of 'exclusion' (in the case when a domain is to be considered as skipped) or 'inclusion', in the case when users wish to include the domain in the solution. Please check the 'as.input' object again.")
        
      } else {
        
        if(length(intersect(x = unique(as.input$cell_type), 
                            y = names(background.networks.list$background.networks))) == 0){
          
          stop("The cell-type names provided in the 'cell_type' column of the 'as.input' object do not match any of the cell-types names provided for 'background.networks' object. Ether revisethe content of 'as.input' data-frame or set to NULL in order to not account for splicing effects.")
          
        } else {
          
          idx2keep <- c()
          for(ii in 1:nrow(as.input)){
            
            ind1 <- which(names(background.networks.list$background.networks)==as.input$cell_type[ii])
            if(length(ind1) == 1){
              
              bn <- background.networks.list$background.networks[[ind1]]
              ind2 <- c(intersect(x = which(bn$pfam_source==as.input$domainID[ii]), 
                                  y = which(bn$gene_source==as.input$proteinID[ii])),
                        intersect(x = which(bn$pfam_target==as.input$domainID[ii]), 
                                  y = which(bn$gene_target==as.input$proteinID[ii])))
              
              if(length(ind2) > 0){idx2keep <- c(idx2keep, ii)}
              
            }
            
          }
          
          if(length(idx2keep) > 0){
            
            print(paste0(length(idx2keep), " domains out of ", nrow(as.input), 
                         " total given in the 'as.input' have been found in the background network."))
            
            as.input <- as.input[idx2keep, ]
            
          } else {
            
            warning("None of the provided domains in the 'as.input' object have been identified in the background network. We will set as.input=NULL and proceed with network optimization without accounting for any splice effects.")
            
            as.input <- NULL
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  
  #### solverPath
  if(!file.exists(solverPath)){
    
    stop("The path to the solver that you provided seem to not exist. Please check your inputs.")
    
  }
  
  
  #### save_res
  if(!is.logical(save_res)){
    
    warning("The 'save_res' object should be logical (TRUE/FALSE). We will set save_res=FALSE.")
    save_res <- FALSE
    
  }
  
  
  #### Lambda's
  if(!is.numeric(lambda1)){
    
    warning("The 'lambda1' parameter should be numeric. We are setting to the default lambda1=10.")
    lambda1 <- 10
    
  }
  
  if(!is.numeric(lambda2)){
    
    warning("The 'lambda2' parameter should be numeric. We are setting to the default lambda2=5.")
    lambda2 <- 5
    
  }
  
  if(!is.numeric(lambda3)){
    
    warning("The 'lambda3' parameter should be numeric. We are setting to the default lambda3=1.")
    lambda3 <- 1
    
  }
  
  if(!is.numeric(lambda4)){
    
    warning("The 'lambda3' parameter should be numeric. We are setting to the default lambda4=1.")
    lambda4 <- 1
    
  }
  
  if(!is.numeric(lambda5)){
    
    warning("The 'lambda5' parameter should be numeric. We are setting to the default lambda5=0.1.")
    lambda5 <- 0.1
    
  }
  
  
  
  #### CPLEX params
  if(!is.numeric(mipgap) ||
     mipgap > 1 ||
     mipgap < 0){
    
    warning("The 'mipgap' parameter should be numeric between 0 and 1. We are stting to the default mipgap=0.05")
    mipgap <- 0.05
    
  }
  
  if(!is.numeric(relgap) ||
     relgap > 1 ||
     relgap < 0){
    
    warning("The 'relgap' parameter should be numeric between 0 and 1. We are stting to the default relgap=0.05.")
    relgap <- 0.05
    
  }
  
  if(!is.numeric(populate) ||
     populate!=round(populate) ||
     populate <= 0){
    
    warning("The 'populate' parameter should be a positive numeric integer. We are stting to the default populate=500.")
    populate <- 500
    
  }
  
  if(!is.numeric(nSolutions) ||
     nSolutions!=round(nSolutions) ||
     nSolutions <= 0){
    
    warning("The 'nSolutions' parameter should be a positive numeric integer. We are stting to the default nSolutions=100.")
    nSolutions <- 100
    
  }
  
  if(!is.numeric(timelimit) ||
     timelimit<=0){
    
    warning("The 'timelimit' parameter should be a positive numeric value. We are stting to the default timelimit=3600 (seconds).")
    timelimit <- 3600
    
  }
  
  if(length(intersect(x = intensity, y = 0:4)) != 1){
    
    warning("The 'intensity' parameter should be numeric value between 0 and 4. We are setting to the default intensity=1.")
    intensity <- 1
    
  }
  
  if(length(intersect(x = replace, y = 0:2)) != 1){
    
    warning("The 'replace' parameter should be numeric value between 0 and 2. We are setting to the default replace=1.")
    replace <- 1
    
  }
  
  if(!is.numeric(threads) ||
     threads!=round(threads) ||
     threads < 0){
    
    warning("The 'threads' parameter should be a positive (or 0) numeric integer. We are setting to the default threads=0.")
    threads <- 0
    
  }
  
  if(!is.numeric(condition) ||
     condition!=round(condition) ||
     condition <= 0){
    
    warning("The 'condition' parameter should be a positive numeric integer. We are setting to the default condition=1.")
    condition <- 1
    
  }
  
  
  
  #### Now do the return object
  all_inputs <- list()
  all_inputs$background.networks.list = background.networks.list
  all_inputs$ccc.input = ccc.input
  all_inputs$tf.input = tf.input
  all_inputs$top.tf = top.tf
  all_inputs$ligand.scores = ligand.scores
  all_inputs$ccc.prob = ccc.prob
  all_inputs$as.input = as.input
  all_inputs$solverPath = solverPath
  all_inputs$lambda1 = lambda1
  all_inputs$lambda2 = lambda2
  all_inputs$lambda3 = lambda3
  all_inputs$lambda4 = lambda4
  all_inputs$lambda5 = lambda5
  all_inputs$mipgap = mipgap
  all_inputs$relgap = relgap
  all_inputs$populate = populate
  all_inputs$nSolutions = nSolutions
  all_inputs$timelimit = timelimit
  all_inputs$intensity = intensity
  all_inputs$replace = replace
  all_inputs$threads = threads
  all_inputs$condition = condition
  all_inputs$save_res = save_res
  
  return(all_inputs)
  
  return(all_inputs)
  
  
}