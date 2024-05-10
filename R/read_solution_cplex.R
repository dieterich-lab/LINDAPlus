read_solution_cplex <- function(variables = variables,
                                background.networks.list = background.networks.list,
                                tf.scores = tf.scores,
                                condition = 1){
  
  ## Network
  cplexSolutionFileName <- paste0("results_", condition, ".txt")
  
  reacVar <- variables$var[which(grepl(pattern = ":reaction ",
                                       x = variables$var_exp))]
  intVar <- variables$var[which(grepl(pattern = ":interaction ",
                                      x = variables$var_exp))]
  ligvar <- variables$var[which(grepl(pattern = "LR:reaction ", 
                                      x = variables$var_exp))]
  varvar <- c(reacVar, intVar, ligvar)
  
  xml_data <- xml2::read_xml(cplexSolutionFileName)
  cplexSolution <- xml2::as_list(xml_data)
  
  #
  cplexSolutionIntAll <- list()
  cplexSolutionReacAll <- list()
  
  sifAll <- list()
  
  for(ii in 2:(length(cplexSolution[[1]]))){
    
    currSolution <- cplexSolution[[1]][[ii]][[4]]
    
    names <- lapply(currSolution, function(x) attr(x, "name"))
    values <- lapply(currSolution, function(x) attr(x, "value"))
    
    mm <- data.frame(name = unlist(names), value = unlist(values), stringsAsFactors = FALSE)
    
    mm <- mm[which(mm$name%in%reacVar), ]
    mm$value <- gsub(pattern = "-", replacement = "", fixed = TRUE, x = mm$value)
    mm$value <- as.numeric(mm$value)
    mm <- mm[which(mm$value>0.99), ]
    
    reacs <- rep("", nrow(mm))
    for(jj in 1:nrow(mm)){
      ind <- which(variables$var==mm$name[jj])
      if(length(ind) > 0){
        reacs[jj] <- variables$var_exp[ind]
      }
    }
    mm$reacs <- reacs
    
    cplexSolutionReacAll[[length(cplexSolutionReacAll)+1]] <- reacs
    
  }
  
  ureacs <- unique(unlist(reacs))
  
  net <- matrix(data = , nrow = length(ureacs), ncol = 5)
  for(ii in 1:length(ureacs)){
    
    cell_type <- strsplit(x = ureacs[ii], split = ":", fixed = TRUE)[[1]][1]
    if(cell_type=="LR"){
      
      ss <- strsplit(x = strsplit(x = ureacs[ii], split = " ", fixed = TRUE)[[1]][2], 
                     split = "=", fixed = TRUE)[[1]][1]
      
      tt <- strsplit(x = strsplit(x = ureacs[ii], split = " ", fixed = TRUE)[[1]][2], 
                     split = "=", fixed = TRUE)[[1]][2]
      
      net[ii, 1] <- cell_type
      net[ii, 2] <- strsplit(x = ss, split = "_", fixed = TRUE)[[1]][2]
      net[ii, 3] <- strsplit(x = tt, split = "_", fixed = TRUE)[[1]][2]
      net[ii, 5] <- paste0(strsplit(x = ss, split = "_", fixed = TRUE)[[1]][1],
                           "=",
                           strsplit(x = tt, split = "_", fixed = TRUE)[[1]][1])
      
      cnt <- 0
      for(jj in 1:length(cplexSolutionReacAll)){
        if(ureacs[ii]%in%cplexSolutionReacAll[[jj]]){cnt <- cnt + 1}
      }
      
      net[ii, 4] <- cnt
      
    } else {
      
      ggs <- strsplit(x = strsplit(x = ureacs[ii], split = " ", fixed = TRUE)[[1]][4], 
                      split = "=", fixed = TRUE)[[1]][1]
      
      ggt <- strsplit(x = strsplit(x = ureacs[ii], split = " ", fixed = TRUE)[[1]][4], 
                      split = "=", fixed = TRUE)[[1]][2]
      
      dds <- strsplit(x = strsplit(x = ureacs[ii], split = " ", fixed = TRUE)[[1]][2], 
                      split = "=", fixed = TRUE)[[1]][1]
      
      ddt <- strsplit(x = strsplit(x = ureacs[ii], split = " ", fixed = TRUE)[[1]][2], 
                      split = "=", fixed = TRUE)[[1]][2]
      
      net[ii, 1] <- cell_type
      net[ii, 2] <- ggs
      net[ii, 3] <- ggt
      net[ii, 5] <- paste0(dds, "=", ddt)
      
      cnt <- 0
      for(jj in 1:length(cplexSolutionReacAll)){
        if(ureacs[ii]%in%cplexSolutionReacAll[[jj]]){cnt <- cnt + 1}
      }
      
      net[ii, 4] <- cnt
      
    }
    
  }
  colnames(net) <- c("Space", "Gene_Source", "Gene_Target", "Weight", "DDI")
  net[which(net[, 1]=="LR"), 1] <- "Extra-Cellular"
  net[, 4] <- as.character(as.numeric(net[, 4])+1)
  net <- as.data.frame(net)
  net$Weight <- as.numeric(net$Weight)
  
  ## Attributes
  allTF <- c()
  for(ii in 1:length(tf.scores)){
    allTF <- c(allTF, tf.scores[[ii]]$tf)
  }
  allTF <- unique(allTF)
  
  uproteins <- unique(c(net$Gene_Source, net$Gene_Target))
  udomains <- c()
  for(ii in 1:nrow(net)){
    udomains <- c(udomains, 
                  c(paste0(net$Gene_Source[ii], "_", strsplit(x = net$DDI[ii], split = "=", fixed = TRUE)[[1]][1]),
                    paste0(net$Gene_Target[ii], "_", strsplit(x = net$DDI[ii], split = "=", fixed = TRUE)[[1]][2]))
    )
  }
  udomains <- unique(udomains)
  
  attr1 <- matrix(data = , nrow = length(uproteins), ncol = 2)
  attr1[, 1] <- uproteins
  attr1[, 2] <- "Protein"
  for(ii in 1:nrow(attr1)){
    
    if(uproteins[ii]%in%background.networks.list$ligands.receptors$ligands){
      attr1[ii, 2] <- "Ligand"
    }
    
    if(uproteins[ii]%in%background.networks.list$ligands.receptors$receptors){
      attr1[ii, 2] <- "Receptor"
    }
    
    if(uproteins[ii]%in%allTF){
      attr1[ii, 2] <- "TF"
    }
    
  }
  
  attr2 <- matrix(data = , nrow = length(udomains), ncol = 2)
  attr2[, 1] <- udomains
  attr2[, 2] <- "Domain"
  
  attributes <- unique(rbind(attr1, attr2))
  colnames(attributes) <- c("node", "attribute")
  attributes <- as.data.frame(attributes)
  
  res <- list()
  res[[length(res)+1]] <- net
  res[[length(res)+1]] <- attributes
  
  names(res) <- c("combined_solutions", "node_attributes")
  
  return(res)
  
  # cplexSolutionFileName <- paste0("results_", condition, ".txt")
  # 
  # xml_data <- xml2::read_xml(cplexSolutionFileName)
  # variables_nodes <- xml2::xml_find_all(xml_data, "//variables/variable")
  # values <- as.numeric(xml2::xml_attr(variables_nodes, "value"))
  # vars <- xml2::xml_attr(variables_nodes, "name")
  # 
  # idx <- which(vars==vars[1])
  # dfList <- list()
  # hop <- idx[2] - idx[1]
  # for(ii in 2:length(idx)){
  #   
  #   df <- matrix(data = , nrow = hop, ncol = 2)
  #   df[, 1] <- vars[(idx[ii]):(idx[ii]+hop-1)]
  #   df[, 2] <- values[(idx[ii]):(idx[ii]+hop-1)]
  #   
  #   dfList[[length(dfList)+1]] <- df
  #   
  # }
  # 
  # 
  # reacVar <- variables$var[which(grepl(pattern = ":reaction ",
  #                                      x = variables$var_exp))]
  # intVar <- variables$var[which(grepl(pattern = ":interaction ",
  #                                     x = variables$var_exp))]
  # ligvar <- variables$var[which(grepl(pattern = "LR:reaction ", 
  #                                     x = variables$var_exp))]
  # 
  # cells <- names(background.networks.list$background.networks)
  # 
  # cplexSolutionData <- XML::xmlParse(cplexSolutionFileName)
  # cplexSolution <- XML::xmlToList(cplexSolutionData)
  # 
  # cplexSolutionIntAll <- list()
  # cplexSolutionReacAll <- list()
  # 
  # sifAll <- list()
  # 
  # for(ii in 2:(length(cplexSolution)-1)){
  #   
  #   currSolution <- cplexSolution[[ii]][[4]]
  #   
  #   mm <- matrix(data = , nrow = length(currSolution), ncol = 3)
  #   for(jj in 1:length(currSolution)){
  #     mm[jj, 1] <- currSolution[[jj]][1]
  #     mm[jj, 2] <- currSolution[[jj]][3]
  #     ind <- which(variables$var==currSolution[[jj]][1])
  #     if(length(ind) > 0){
  #       mm[jj, 3] <- variables$var_exp[which(variables$var==currSolution[[jj]][1])]
  #     }
  #   }
  #   mm <- mm[complete.cases(mm), ]
  #   
  #   varvar <- unlist(lapply(currSolution, '[', 1))
  #   valval <- as.numeric(unlist(lapply(currSolution, '[', 3)))
  #   
  #   idxReac <- intersect(x = which(varvar%in%reacVar), y = which(valval>=0.99))
  #   if(length(idxReac)>0){
  #     
  #     temp <- list()
  #     
  #     for(kk in 1:length(cells)){
  #       
  #       reactions <-
  #         sapply(strsplit(x =
  #                           variables$var_exp[
  #                             intersect(x = which(variables$var%in%varvar[idxReac]), 
  #                                       y = which(grepl(pattern = paste0(cells[kk], ":"), x = variables$var_exp)))],
  #                         split = " ", fixed = TRUE), '[', 2)
  #       interactions <-
  #         sapply(strsplit(x =
  #                           variables$var_exp[
  #                             intersect(x = which(variables$var%in%varvar[idxReac]), 
  #                                       y = which(grepl(pattern = paste0(cells[kk], ":"), x = variables$var_exp)))],
  #                         split = " ", fixed = TRUE), '[', 4)
  #       
  #       # uInt <- unique(interactions)
  #       uInt <- interactions
  #       currSIF <- matrix(data = , nrow = length(uInt), ncol = 4)
  #       currSIF[, 1] <- sapply(strsplit(x = uInt,
  #                                       split = "=",
  #                                       fixed = TRUE),
  #                              '[',
  #                              1)
  #       currSIF[, 2] <- "1"
  #       currSIF[, 3] <- sapply(strsplit(x = uInt,
  #                                       split = "=",
  #                                       fixed = TRUE),
  #                              '[',
  #                              2)
  #       currSIF[, 4] <- reactions
  #       # for(jj in seq_len(length(uInt))){
  #       #   
  #       #   idx <-
  #       #     which(grepl(pattern =
  #       #                   paste0(" ",
  #       #                          uInt[jj]),
  #       #                 x = variables$var_exp[
  #       #                   which(variables$var%in%varvar[idxReac])],
  #       #                 fixed = TRUE))
  #       #   currSIF[jj, 4] <- paste0(reactions[idx], collapse = "; ")
  #       #   
  #       # }
  #       
  #       temp[[length(temp)+1]] <- currSIF
  #       
  #     }
  #     
  #     idxLig <- intersect(x = which(varvar%in%ligvar), y = which(valval>=0.99))
  #     lig_reactions <- variables$var_exp[which(variables$var%in%varvar[idxLig])]
  #     lig_reactions <- sapply(strsplit(x = lig_reactions, split = " ", fixed = TRUE), "[", 2)
  #     lig_mat <- matrix(data = , nrow = length(lig_reactions), ncol = 4)
  #     ss <- sapply(strsplit(x = lig_reactions, split = "=", fixed = TRUE), "[", 1)
  #     tt <- sapply(strsplit(x = lig_reactions, split = "=", fixed = TRUE), "[", 2)
  #     lig_mat[, 1] <- sapply(strsplit(x = ss, split = "_", fixed = TRUE), "[", 2)
  #     lig_mat[, 2] <- "1"
  #     lig_mat[, 3] <- sapply(strsplit(x = tt, split = "_", fixed = TRUE), "[", 2)
  #     lig_mat[, 4] <- paste0(sapply(strsplit(x = ss, split = "_", fixed = TRUE), "[", 1),
  #                            "=",
  #                            sapply(strsplit(x = tt, split = "_", fixed = TRUE), "[", 1))
  #     
  #     temp[[length(temp)+1]] <- lig_mat
  #     names(temp) <- c(cells, "ligand_receptors")
  #     
  #     sifAll[[length(sifAll)+1]] <- temp
  #     
  #   }
  #   
  # }
  # names(sifAll) <- paste0("Solution-", 1:length(sifAll))
  # separate_solutions <- sifAll
  # 
  # cases <- length(sifAll$`Solution-1`)
  # combined_solutions <- vector("list", length = length(cases))
  # for(kk in 1:cases){
  #   
  #   for(ii in seq_len(length(sifAll))){
  #     
  #     if(ii==1){
  #       
  #       combSIF <- sifAll[[ii]][[kk]][, 1:3]
  #       if(class(combSIF)[1]=="character"){
  #         combSIF <- t(as.matrix(combSIF))
  #       }
  #       
  #     } else {
  #       
  #       currSIF <- sifAll[[ii]][[kk]][, 1:3]
  #       if(class(currSIF)[1]=="character"){
  #         currSIF <- t(as.matrix(currSIF))
  #       }
  #       for(jj in seq_len(nrow(currSIF))){
  #         
  #         idx1 <- which(combSIF[, 1]==currSIF[jj, 1])
  #         idx2 <- which(combSIF[, 3]==currSIF[jj, 3])
  #         idx <- intersect(x = idx1, y = idx2)
  #         if(length(idx)>0){
  #           
  #           for(zz in 1:length(idx)){combSIF[idx[zz], 2] <- as.character(as.numeric(combSIF[idx[zz], 2])+1)}
  #           
  #         } else {
  #           
  #           combSIF <- rbind(combSIF, t(as.matrix(currSIF[jj, ])))
  #           
  #         }
  #         
  #       }
  #       
  #     }
  #     
  #   }
  #   
  #   domains <- rep("", nrow(combSIF))
  #   for(ii in seq_len(nrow(combSIF))){
  #     
  #     ss <- combSIF[ii, 1]
  #     tt <- combSIF[ii, 3]
  #     
  #     ud <- c()
  #     for(jj in seq_len(length(sifAll))){
  #       
  #       idx1 <- which(sifAll[[jj]][[kk]][, 1]==ss)
  #       idx2 <- which(sifAll[[jj]][[kk]][, 3]==tt)
  #       idx <- intersect(x = idx1, y = idx2)
  #       if(length(idx)>0){
  #         
  #         ud <- c(ud, strsplit(x = sifAll[[jj]][[kk]][idx, 4],
  #                              split = "; ",
  #                              fixed = TRUE)[[1]])
  #         
  #       }
  #       
  #     }
  #     
  #     domains[ii] <- paste0(unique(ud), collapse = "; ")
  #   }
  #   
  #   sif <- matrix(data = , nrow = nrow(combSIF), ncol = 4)
  #   sif[, 1:3] <- combSIF
  #   sif[, 4] <- domains
  #   
  #   if(nrow(sif) > 0){
  #     colnames(sif) <- c("source", "weight", "target", "reaction")
  #   }
  #   
  #   combined_solutions[[kk]] <- sif
  #   
  # }
  # names(combined_solutions) <- c(cells, "ligand_receptors")
  # 
  # attributes <- list()
  # for(ii in 1:length(background.networks.list$background.networks)){
  #   
  #   nodes <- unique(c(background.networks.list$background.networks[[ii]]$pfam_source,
  #                     background.networks.list$background.networks[[ii]]$pfam_target,
  #                     background.networks.list$background.networks[[ii]]$gene_source,
  #                     background.networks.list$background.networks[[ii]]$gene_target))
  #   
  #   mm <- matrix(data = , nrow = length(nodes), ncol = 2)
  #   mm[, 1] <- nodes
  #   mm[, 2] <- "protein"
  #   
  #   ind <- which(nodes%in%c(background.networks.list$background.networks[[ii]]$pfam_source,
  #                           background.networks.list$background.networks[[ii]]$pfam_target))
  #   if(length(ind) > 0){mm[ind, 2] <- "domain"}
  #   
  #   ind <- which(nodes%in%tf.scores[[ii]]$tf)
  #   if(length(ind) > 0){mm[ind, 2] <- "tf"}
  #   
  #   ind <- which(nodes%in%background.networks.list$ligands.receptors$ligands)
  #   if(length(ind) > 0){mm[ind, 2] <- "ligand"}
  #   
  #   ind <- which(nodes%in%background.networks.list$ligands.receptors$receptors)
  #   if(length(ind) > 0){mm[ind, 2] <- "receptor"}
  #   
  #   mm <- mm[c(which(mm[, 2]=="receptor"), which(mm[, 2]=="protein"), 
  #              which(mm[, 2]=="tf"), which(mm[, 2]=="ligand"), which(mm[, 2]=="domain")), ]
  #   
  #   colnames(mm) <- c("node", "attribute")
  #   
  #   attributes[[length(attributes)+1]] <- mm
  #   
  # }
  # 
  # lr <- combined_solutions$ligand_receptors
  # nodes <- unique(c(lr[, 1], lr[, 3]))
  # mm <- matrix(data = , nrow = length(nodes), ncol = 2)
  # mm[, 1] <- nodes
  # for(ii in 1:length(nodes)){
  #   if(nodes[ii] %in% lr[, 1]){
  #     mm[ii, 2] <- "ligand"
  #   } else {
  #     mm[ii, 2] <- "receptor"
  #   }
  # }
  # colnames(mm) <- c("node", "attribute")
  # attributes[[length(attributes)+1]] <- mm
  # 
  # names(attributes) <- c(names(background.networks.list$background.networks), "ligand-receptors")
  # 
  # res <- list()
  # res[[length(res)+1]] <- combined_solutions
  # res[[length(res)+1]] <- separate_solutions
  # res[[length(res)+1]] <- attributes
  # names(res) <- c("combined_solutions", "separate_solutions", "node_attributes")
  # 
  # return(res)
  
}
