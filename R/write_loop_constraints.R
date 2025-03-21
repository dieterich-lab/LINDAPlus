write_loop_constraints <- function(variables = variables,
                                   background.networks.list = background.networks.list,
                                   constraits.parallel.writing = constraits.parallel.writing) {
  
  ligands <- background.networks.list$ligands.receptors$ligands
  receptors <- background.networks.list$ligands.receptors$receptors
  cells <- names(background.networks.list$background.networks)
  
  lr <- variables$var_exp[grepl("LR:interaction ", variables$var_exp, fixed = TRUE)]
  lr <- sapply(strsplit(lr, " "), function(x) x[2])
  
  process_cell_constraints <- function(ii) {
    background.network <- background.networks.list$background.networks[[ii]]
    sif <- unique(background.network[, c("gene_source", "gene_target")])
    lr_mask <- paste0(sif[, 1], "=", sif[, 2]) %in% lr
    sif <- sif[!lr_mask, ]
    
    speciesVar <- variables$var[grepl(paste0(cells[ii], ":node "), variables$var_exp, fixed = TRUE)]
    speciesExp <- variables$var_exp[grepl(paste0(cells[ii], ":node "), variables$var_exp, fixed = TRUE)]
    
    reacVar <- variables$var[grepl(paste0(cells[ii], ":interaction "), variables$var_exp, fixed = TRUE)]
    distVar <- variables$var[grepl(paste0(cells[ii], ":dist "), variables$var_exp, fixed = TRUE)]
    
    cc1 <- paste0(speciesVar, " - ", distVar, " <= 0")
    cc2 <- paste0(distVar, " <= 10001")
    
    cc3 <- vector("character", nrow(sif))
    for (jj in seq_len(nrow(sif))) {
      var1_idx <- match(paste0(cells[ii], ":dist ", sif[jj, 2]), variables$var_exp)
      var2_idx <- match(paste0(cells[ii], ":dist ", sif[jj, 1]), variables$var_exp)
      var3_idx <- match(paste0(cells[ii], ":interaction ", sif[jj, 1], "=", sif[jj, 2]), variables$var_exp)
      
      cc3[jj] <- paste0(variables$var[var1_idx],
                        " - ",
                        variables$var[var2_idx],
                        " + 10001 ",
                        variables$var[var3_idx],
                        " <= ",
                        -(1 - 10001))
    }
    
    return(unique(c(cc1, cc2, cc3)))
  }
  
  if (constraits.parallel.writing) {
    constraints_list <- mclapply(seq_along(cells), process_cell_constraints, mc.cores = length(cells))
  } else {
    constraints_list <- lapply(seq_along(cells), process_cell_constraints)
  }
  
  return(unlist(constraints_list))
}
