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
  
  # xml_data <- xml2::read_xml(cplexSolutionFileName)
  # cplexSolution <- xml2::as_list(xml_data)
  
  process_solution <- function(solution) {
    currSolution <- xml2::xml_find_all(solution, ".//variable")
    names <- sapply(currSolution, function(x) xml2::xml_attr(x, "name"))
    values <- sapply(currSolution, function(x) xml2::xml_attr(x, "value"))
    list(names = names, values = values)
  }
  
  xml_data <- xml2::read_xml(cplexSolutionFileName, options = "COMPACT")
  solutions <- xml2::xml_find_all(xml_data, ".//CPLEXSolution")
  results <- lapply(solutions[2:length(solutions)], process_solution)
  # save(results, file = paste0("result_", condition, ".RData"))
  
  #
  cplexSolutionIntAll <- list()
  cplexSolutionReacAll <- list()
  
  sifAll <- list()
  
  for(ii in 1:(length(results))){
    
    currSolution <- results[[ii]]
    
    names <- currSolution$names
    values <- currSolution$values
    
    MM <- data.frame(name = names, value = values, stringsAsFactors = FALSE)
    
    mm <- MM[which(MM$name%in%reacVar), ]
    mm$value <- gsub(pattern = "-", replacement = "", fixed = TRUE, x = mm$value)
    mm$value <- as.numeric(mm$value)
    mm <- mm[intersect(x = which(mm$value>0.99), y = which(mm$value<1.01)), ]
    
    reacs <- rep("", nrow(mm))
    for(jj in 1:nrow(mm)){
      ind <- which(variables$var==mm$name[jj])
      if(length(ind) > 0){
        reacs[jj] <- variables$var_exp[ind]
      }
    }
    mm$reacs <- reacs
    
    mm <- MM[which(MM$name%in%intVar), ]
    mm$value <- gsub(pattern = "-", replacement = "", fixed = TRUE, x = mm$value)
    mm$value <- as.numeric(mm$value)
    mm <- mm[intersect(x = which(mm$value>0.99), y = which(mm$value<1.01)), ]
    
    intPPI <- variables$var_exp[which(variables$var%in%mm$name)]
    intPPI <- paste0(sapply(strsplit(x = intPPI, split = ":", fixed = TRUE), "[", 1), ":", 
                     sapply(strsplit(x = intPPI, split = " ", fixed = TRUE), "[", 2))
    
    intDDI <- paste0(sapply(strsplit(x = reacs, split = ":", fixed = TRUE), "[", 1), ":", 
                     sapply(strsplit(x = reacs, split = " ", fixed = TRUE), "[", 4))
    indNA <- which(intDDI %in% paste0(c("LR", names(background.networks.list$background.networks)), ":NA"))
    for(zz in 1:length(indNA)){
      
      curr <- reacs[indNA[zz]]
      curr <- strsplit(x = curr, split = " ", fixed = TRUE)[[1]][2]
      currSS <- strsplit(x = strsplit(x = curr, split = "=", fixed = TRUE)[[1]][1], 
                         split = "_", fixed = TRUE)[[1]][2]
      currTT <- strsplit(x = strsplit(x = curr, split = "=", fixed = TRUE)[[1]][2], 
                         split = "_", fixed = TRUE)[[1]][2]
      intDDI[indNA[zz]] <- paste0(strsplit(reacs[indNA[zz]], split = ":", fixed = TRUE)[[1]][1], 
                                  ":", currSS, "=", currTT)
      
    }
    idx2keep <- which(intDDI %in% intPPI)
    
    cplexSolutionReacAll[[length(cplexSolutionReacAll)+1]] <- reacs[idx2keep]
    
  }
  
  ureacs <- unique(unlist(cplexSolutionReacAll))
  
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
  ct <- setdiff(x = unique(net$Space), y = "Extra-Cellular")
  ct <- ct[order(ct)]
  ind2order <- which(net$Space=="Extra-Cellular")
  for(ii in 1:length(ct)){
    ind2order <- c(ind2order, which(net$Space==ct[ii]))
  }
  net <- net[ind2order, ]
  ind2rem <- which(net$Gene_Source%in%c("PSEUDOLIGAND", "PSEUDORECEPTOR", "PSEUDOPROTEIN", "PSEUDOTF"))
  if(length(ind2rem) > 0){net <- net[-ind2rem, ]}
  
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
  
}
