write_constraints_2 <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  constraints <- c()
  
  cell_types <- names(background.networks.list$background.networks)
  
  for(ii in 1:length(cell_types)){
    
    var <- variables$var[which(grepl(pattern = paste0(cell_types[ii], ":"), x = variables$var_exp, fixed = TRUE))]
    var_exp <- variables$var_exp[which(grepl(pattern = paste0(cell_types[ii], ":"), x = variables$var_exp, fixed = TRUE))]
    var_exp <- gsub(pattern = paste0(cell_types[ii], ":"), replacement = "", x = var_exp, fixed = TRUE)
    
    # Interactions
    cc1 <- c()
    idx <- which(grepl(pattern = "interaction ", x = var_exp, fixed = TRUE))
    ss <- sapply(strsplit(x = gsub(pattern = "interaction ", replacement = "", x = var_exp[idx], fixed = TRUE), 
                          split = "=", fixed = TRUE), "[", 1)
    tt <- sapply(strsplit(x = gsub(pattern = "interaction ", replacement = "", x = var_exp[idx], fixed = TRUE), 
                          split = "=", fixed = TRUE), "[", 2)
    for(jj in 1:length(idx)){
      
      cc1 <- c(cc1, paste0(var[which(var_exp==paste0("node ", ss[jj]))], 
                           " + ",
                           var[which(var_exp==paste0("node ", tt[jj]))],
                           " - 2 ", var[idx[jj]],
                           " >= 0"))
    }
    
    # Reactions
    cc2 <- c()
    idx <- which(grepl(pattern = "reaction ", x = var_exp, fixed = TRUE))
    ddi <- sapply(strsplit(x = var_exp[idx], split = " ", fixed = TRUE), "[", 2)
    ppi <- sapply(strsplit(x = var_exp[idx], split = " ", fixed = TRUE), "[", 4)
    dd1 <- sapply(strsplit(x = ddi, split = "=", fixed = TRUE), "[", 1)
    dd2 <- sapply(strsplit(x = ddi, split = "=", fixed = TRUE), "[", 2)
    pp1 <- sapply(strsplit(x = ppi, split = "=", fixed = TRUE), "[", 1)
    pp2 <- sapply(strsplit(x = ppi, split = "=", fixed = TRUE), "[", 2)
    for(jj in 1:length(idx)){
      
      ind1 <- which(var_exp==paste0("domain ", dd1[jj], " of protein ", pp1[jj]))
      ind2 <- which(var_exp==paste0("domain ", dd2[jj], " of protein ", pp2[jj]))
      
      cc2 <- c(cc2, paste0(var[ind1], " + ", var[ind2], " - 2 ", var[idx[jj]], " >= 0"))
      
      
    }
    
    constraints <- c(constraints, cc1, cc2)
    
  }
  
  return(constraints)
  
}