checkInputs <- function(background.networks.list = background.networks.list,
                        ligand.scores = ligand.scores,
                        tf.scores = tf.scores,
                        solverPath = solverPath,
                        top = top,
                        alpha = alpha,
                        beta = beta,
                        gamma = gamma,
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
  
  # background.networks.list
  if((class(background.networks.list) != "list") ||
     (length(background.networks.list)!= 2) || 
     (length(intersect(x = names(background.networks.list), y = c("background.networks", "ligands.receptors")))<2)){
    
    stop("Error: Object 'background.networks.list' should be a list with a length of 2 named 'background.networks' and 'ligands.receptors'. Please check your inputs.")
    
  } else {
    
    obj1 <- background.networks.list$background.networks
    obj2 <- background.networks.list$ligands.receptors
    
    #obj1
    cc <- c()
    for(ii in 1:length(obj1)){cc <- c(cc, class(obj1[[ii]]))}
    if((is.null(names(obj1))) ||
       all(cc != "data.frame")){
      stop("Error: The 'background.networks' should be named data-frames. Please check your inputs.")
    } else {
      
      for(ii in 1:length(obj1)){
        
        df <- obj1[[ii]]
        if(length(intersect(x = colnames(df), y = c("pfam_source", "pfam_target", "gene_source", "gene_target")))<4){
          
          stop("Error: The 'background.network' object should be a data-frame with at least 4
              columns and have at lease c('pfam_source', 'pfam_target', 
              'gene_source', 'gene_target') as column ID's. Please check your inputs.")
        }
        
      }
    }
    
    #obj2
    cc <- c()
    for(ii in 1:length(obj2)){cc <- c(cc, class(obj2[[ii]]))}
    if(length(intersect(x = names(obj2), y = c("ligands", "receptors")))<2 ||
       all(cc != "character")){
      stop("Error: The 'ligands.receptors' should be a list of characters with a length of 2 named 'ligands' and 'receptors'. Please check your inputs.")
    }
    
  }
  
  
  # tf.scores
  if((class(tf.scores) != "list") ||
     (length(tf.scores)!= length(background.networks.list$background.networks)) || 
     (length(intersect(x = names(background.networks.list$background.networks), y = names(tf.scores)))!=length(tf.scores))){
    
    stop("Error: The 'tf.scores' object should be a list of data-frames with the 
         same length and named the same as the 'background.networks.list$background.networks' 
         object. Please check your inputs.")
    
  } else {
    
    cc <- c()
    for(ii in 1:length(tf.scores)){cc <- c(cc, class(tf.scores[[ii]]))}
    if(all(cc != "data.frame")){
      stop("Error: The 'tf.scores' object should be a list of data-frames with the 
           same length and named the same as the 'background.networks.list$background.networks' 
           object. Please check your inputs.")
    }
    
  }
  
  if((!is.numeric(top)) ||
     (length(intersect(x = names(top), y = names(background.networks.list$background.networks)))<length(top))){
    
    stop("Error: The 'top' parameter should be a named numeric vector with cell-types given as names. Please check your inputs.")
    
  } else {
    
    temp <- list()
    for(ii in 1:length(top)){
      
      curr <- tf.scores[[which(names(tf.scores)==names(top)[ii])]]
      
      if(top[ii] > nrow(curr)){
        
        warning(paste0("WARNING: There a re more given top TF's for cell-type ",
                       names(top)[ii], " than there are TF's given. All the given
                       TF's will be considered as significant."))
        
        curr$score <- 1
        temp[[length(temp)+1]] <- curr
        
      } else {
        
        curr <- curr[order(abs(as.numeric(curr$score)), decreasing = TRUE), ]
        curr$score <- 0
        curr$score[1:top[ii]] <- 1
        temp[[length(temp)+1]] <- curr
        
      }
      
    }
    names(temp) <- names(top)
    tf.scores <- temp
    
  }
  
  if(!is.numeric(alpha)){
    
    warning("The 'alpha' parameter should be numeric. We are setting to the default alpha=10.")
    alpha <- 10
    
  }
  
  if(!is.numeric(beta)){
    
    warning("The 'beta' parameter should be numeric. We are setting to the default beta=5.")
    beta <- 5
    
  }
  
  if(!is.numeric(gamma)){
    
    warning("The 'gamma' parameter should be numeric. We are setting to the default beta=0.1.")
    gamma <- 0.1
    
  }
  
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
    
    warning("The 'nSolutions' parameter should be a positive numeric integer.
         We are stting to the default nSolutions=100.")
    nSolutions <- 100
    
  }
  
  if(!is.numeric(timelimit) ||
     timelimit<=0){
    
    warning("The 'timelimit' parameter should be a positive numeric value. We are stting to the default timelimit=3600 (seconds).")
    timelimit <- 3600
    
  }
  
  if(length(intersect(x = intensity, y = 0:4)) != 1){
    
    warning("The 'intensity' parameter should be numeric value between 0 and 4.
             We are setting to the default intensity=1.")
    intensity <- 1
    
  }
  
  if(length(intersect(x = replace, y = 0:2)) != 1){
    
    warning("The 'replace' parameter should be numeric value between 0 and 2.
             We are setting to the default replace=1.")
    replace <- 1
    
  }
  
  if(!is.numeric(threads) ||
     threads!=round(threads) ||
     threads < 0){
    
    warning("The 'threads' parameter should be a positive (or 0) numeric integer. 
             We are setting to the default threads=0.")
    threads <- 0
    
  }
  
  if(!is.numeric(condition) ||
     condition!=round(condition) ||
     condition <= 0){
    
    warning("The 'condition' parameter should be a positive numeric integer.
             We are setting to the default condition=1.")
    condition <- 1
    
  }
  
  # solverPath
  if(!file.exists(solverPath)){
    
    stop("The path to the solver that you provided seem to not exist. Please check your inputs.")
    
  }
  
  # save_res
  if(!is.logical(save_res)){
    
    warning("The 'save_res' object should be logical (TRUE/FALSE). We will set save_res=FALSE.")
    save_res <- FALSE
    
  }
  
  all_inputs <- list()
  all_inputs$background.networks.list = background.networks.list
  all_inputs$ligand.scores = ligand.scores
  all_inputs$tf.scores = tf.scores
  all_inputs$solverPath = solverPath
  all_inputs$top = top
  all_inputs$alpha = alpha
  all_inputs$beta = beta
  all_inputs$gamma = gamma
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
  
}

