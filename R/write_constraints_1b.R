write_constraints_1b <- function(variables = variables,
                                background.networks.list = background.networks.list){
  
  constraints <- c()
  
  receptors <- background.networks.list$ligands.receptors$receptors
  ligands <- background.networks.list$ligands.receptors$ligands
  
  background.network.list <- background.networks.list$background.networks
  
  # There should be at least one Receptor activated on each cell
  cell_types <- names(background.networks.list$background.networks)
  
  # Make the connection between LR: and CellType: interactions
  # If there is one LR interaction active, there should be at least one
  # CellType interaction active
  lr_int <- variables$var_exp[which(grepl(pattern = "LR:interaction", x = variables$var_exp, fixed = TRUE))]
  for(ii in 1:length(lr_int)){
    
    curr <- lr_int[ii]
    vv1 <- which(variables$var_exp==curr)
    vv <- c()
    for(jj in 1:length(cell_types)){
      vv <- c(vv, which(variables$var_exp==gsub(pattern = "LR:", replacement = paste0(cell_types[jj], ":"), x = curr, fixed = TRUE)))
    }
    
    con1 <- paste0(variables$var[vv1], " - ", paste0(variables$var[vv], collapse = " - "), " <= 0")
    con2 <- paste0(variables$var[vv1], " - ", variables$var[vv], " >= 0")
    
    constraints <- c(constraints, c(con1, con2))
    
  }
  
  # Make the connection between LR: and CellType: interactions
  # If there is one protein active, there should be at least one
  # CellType protein active
  lr_prot <- c(variables$var_exp[which(grepl(pattern = "LR:receptor", x = variables$var_exp, fixed = TRUE))],
               variables$var_exp[which(grepl(pattern = "LR:ligand", x = variables$var_exp, fixed = TRUE))])
  unodes <- unique(sapply(strsplit(x = lr_prot, split = " ", fixed = TRUE), "[", 2))
  for(ii in 1:length(lr_prot)){
    
    curr <- lr_prot[ii]
    vv1 <- which(variables$var_exp==curr)
    vv <- c()
    for(jj in 1:length(cell_types)){
      if(grepl(pattern = "receptor", x = curr, fixed = TRUE)){
        vv <- c(vv, which(variables$var_exp==gsub(pattern = "LR:receptor", replacement = paste0(cell_types[jj], ":node"), x = curr, fixed = TRUE)))
      } else {
        vv <- c(vv, which(variables$var_exp==gsub(pattern = "LR:ligand", replacement = paste0(cell_types[jj], ":node"), x = curr, fixed = TRUE)))
      }
    }
    
    con1 <- paste0(variables$var[vv1], " - ", paste0(variables$var[vv], collapse = " - "), " <= 0")
    con2 <- paste0(variables$var[vv1], " - ", variables$var[vv], " >= 0")
    
    constraints <- c(constraints, c(con1, con2))
    
  }
  
  # Make the connection between LR: and CellType: interactions
  # If there is one domain active, there should be at least one
  # CellType domain active
  lr_prot <- c(variables$var_exp[which(grepl(pattern = "LR:domain_receptor", x = variables$var_exp, fixed = TRUE))],
               variables$var_exp[which(grepl(pattern = "LR:domain_ligand", x = variables$var_exp, fixed = TRUE))])
  unodes <- unique(sapply(strsplit(x = lr_prot, split = " ", fixed = TRUE), "[", 2))
  for(ii in 1:length(lr_prot)){
    
    curr <- lr_prot[ii]
    vv1 <- which(variables$var_exp==curr)
    vv <- c()
    dd <- strsplit(x = strsplit(x = curr, split = " ", fixed = TRUE)[[1]][2], split = "_", fixed = TRUE)[[1]][1]
    pp <- strsplit(x = strsplit(x = curr, split = " ", fixed = TRUE)[[1]][2], split = "_", fixed = TRUE)[[1]][2]
    for(jj in 1:length(cell_types)){
      vv <- c(vv, which(variables$var_exp==paste0(cell_types[jj], ":domain ", dd, " of protein ", pp)))
    }
    
    con1 <- paste0(variables$var[vv1], " - ", paste0(variables$var[vv], collapse = " - "), " <= 0")
    con2 <- paste0(variables$var[vv1], " - ", variables$var[vv], " >= 0")
    
    constraints <- c(constraints, c(con1, con2))
    
  }
  
  
  # A receptor in cell-type k is present in the solution, then
  # there should be at least one functional interaction with a
  # ligand upstream of it

  # Protein Level
  for(ii in 1:length(cell_types)){

    var <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
    var_exp <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]

    bg <- background.networks.list$background.networks[[ii]]
    

    lr1 <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR")]
    lr_exp1 <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR")]
    lr2 <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
    lr_exp2 <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
    ind <- which(gsub(pattern = cell_types[ii], replacement = "LR", x = lr_exp2) %in% lr_exp1)
    lr <- lr2[ind]
    lr_exp <- lr_exp2[ind]
    for(jj in 1:length(receptors)){

      idx <- intersect(x = which(bg$gene_source%in%ligands), y = which(bg$gene_target%in%receptors[jj]))
      if(length(idx) > 0){

        tmp <- bg[idx, ]
        ind1 <- which(var_exp==paste0(cell_types[ii], ":node ", receptors[jj]))

        ind2 <- which(lr_exp%in%paste0(cell_types[ii], ":interaction ", tmp$gene_source, "=", tmp$gene_target))

        cc2 <- paste0(lr[ind2], collapse = " - ")
        cc2 <- paste0(var[ind1], " - ", cc2, " <= 0")
        constraints <- c(constraints, cc2)

        cc22 <- paste0(var[ind1], " - ", lr[ind2], " >= 0")
        constraints <- c(constraints, cc22)

      }

    }

  }

  # Domain Level
  for(ii in 1:length(cell_types)){

    var <- variables$var[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]
    var_exp <- variables$var_exp[which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii])]

    bg <- background.networks.list$background.networks[[ii]]

    lr1 <- variables$var[intersect(x = which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR"), 
                                   y = which(grepl(pattern = "reaction PSEUDODOMAINLR", x = variables$var_exp, fixed = TRUE)))]
    lr_exp1 <- variables$var_exp[intersect(x = which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)=="LR"), 
                                           y = which(grepl(pattern = "reaction PSEUDODOMAINLR", x = variables$var_exp, fixed = TRUE)))]
    lr2 <- variables$var[intersect(x = which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii]), 
                                   y = which(grepl(pattern = "reaction PSEUDODOMAINLR", x = variables$var_exp, fixed = TRUE)))]
    lr_exp2 <- variables$var_exp[intersect(x = which(sapply(strsplit(x = variables$var_exp, split = ":", fixed = TRUE), "[", 1)==cell_types[ii]), 
                                           y = which(grepl(pattern = "reaction PSEUDODOMAINLR", x = variables$var_exp, fixed = TRUE)))]
    for(jj in 1:length(receptors)){

      idx <- intersect(x = which(bg$gene_source%in%ligands), y = which(bg$gene_target%in%receptors[jj]))
      if(length(idx) > 0){

        tmp <- bg[idx, ]
        udomains <- unique(paste0(tmp$pfam_target))

        for(kk in 1:length(udomains)){

          tmptmp <- tmp[which(tmp$pfam_target==udomains[kk]), ]

          ind1 <- which(lr_exp2%in%paste0(cell_types[ii], ":reaction ", tmptmp$pfam_source, "=", tmptmp$pfam_target,
                                          " of ", tmptmp$gene_source, "=", tmptmp$gene_target))

          ind2 <- which(var_exp==paste0(cell_types[ii], ":domain ", udomains[kk],
                                        " of protein ", receptors[jj]))


          cc3 <- paste0(lr2[ind1], collapse = " - ")
          cc3 <- paste0(var[ind2], " - ", cc3, " <= 0")
          constraints <- c(constraints, cc3)

          cc33 <- paste0(var[ind2], " - ", lr2[ind1], " >= 0")
          constraints <- c(constraints, cc33)

        }

      }

    }

  }
  
  ind2rem <- which(grepl(pattern = " -  <= 0", x = constraints, fixed = TRUE))
  if(length(ind2rem) > 0){constraints <- constraints[-ind2rem]}
  
  ind2rem <- which(grepl(pattern = " -  >= 0", x = constraints, fixed = TRUE))
  if(length(ind2rem) > 0){constraints <- constraints[-ind2rem]}
  
  return(constraints)
  
}