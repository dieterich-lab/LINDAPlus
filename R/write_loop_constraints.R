write_loop_constraints <- function(variables = variables,
                                   background.networks.list = background.networks.list){
  
  ligands <- background.networks.list$ligands.receptors$ligands
  receptors <- background.networks.list$ligands.receptors$receptors
  cells <- names(background.networks.list$background.networks)
  
  lr <- variables$var_exp[which(grepl(pattern = "LR:interaction ", 
                                      x = variables$var_exp, 
                                      fixed = TRUE))]
  lr <- sapply(strsplit(x = lr, split = " ", fixed = TRUE), "[", 2)
  
  constraints <- c()
  for(ii in 1:length(background.networks.list$background.networks)){
    
    background.network <- background.networks.list$background.networks[[ii]]
    
    sif <- unique(background.network[, c("gene_source", "gene_target")])
    
    idx2rem <- which(paste0(sif[, 1], "=", sif[, 2])%in%lr)
    if(length(idx2rem) > 0){
      sif <- sif[-idx2rem, ]
    }
    
    species <- unique(c(sif[, 1], sif[, 2]))
    
    speciesVar <- variables$var[grepl(pattern = paste0(cells[ii], ":node "),
                                      x = variables$var_exp, 
                                      fixed = TRUE)]
    speciesExp <- variables$var_exp[grepl(pattern = paste0(cells[ii], ":node "),
                                          x = variables$var_exp, 
                                          fixed = TRUE)]
    reacVar <- variables$var[grepl(pattern = paste0(cells[ii], ":interaction "),
                                   x = variables$var_exp, 
                                   fixed = TRUE)]
    distVar <- variables$var[grepl(pattern = paste0(cells[ii], ":dist "),
                                   x = variables$var_exp, 
                                   fixed = TRUE)]
    
    cc1 <- paste0(speciesVar, " - ", distVar, " <= 0")
    
    cc2 <- paste0(distVar, " <= ", 10001)
    
    cc3 <- rep("", nrow(sif))
    for(jj in seq_len(nrow(sif))){
      
      var1 <- variables$var[which(variables$var_exp==paste0(cells[ii], ":dist ",
                                                            sif[jj, 2]))]
      var2 <- variables$var[which(variables$var_exp==paste0(cells[ii], ":dist ",
                                                            sif[jj, 1]))]
      var3 <- variables$var[which(variables$var_exp==paste0(cells[ii], ":interaction ",
                                                            sif[jj, 1],
                                                            "=",
                                                            sif[jj, 2]))]
      
      cc3[jj] <- paste0(var1,
                        " - ",
                        var2,
                        " + 10001 ",
                        var3,
                        " <= ",
                        -(1-10001))
      
    }
    
    constraints <- unique(c(constraints, c(cc1, cc2, cc3)))
    
  }
  
  
  
  return(constraints)
  
}
