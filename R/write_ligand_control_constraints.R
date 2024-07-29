write_ligand_control_constraints <-  function(variables = variables,
                                              background.networks.list = background.networks.list){
  
  ind <- which(grepl(pattern = "reaction PSEUDODOMAINTF=PSEUDODOMAINTF", 
                     x = variables$var_exp, 
                     fixed = TRUE))
  
  uligands <- variables$var_exp[which(grepl(pattern = "reaction PSEUDODOMAINTF=PSEUDODOMAINTF", 
                                            x = variables$var_exp, 
                                            fixed = TRUE))]
  
  get_last_element <- function(x) {
    split_elements <- strsplit(x, " ")[[1]]
    return(tail(split_elements, 1))
  }
  uligands <- as.character(sapply(uligands, get_last_element))
  
  get_last_element2 <- function(x) {
    split_elements <- strsplit(x, "=")[[1]]
    return(tail(split_elements, 1))
  }
  uligands_ctrl <- as.character(sapply(uligands, get_last_element2))
  uligands <- unique(as.character(sapply(uligands, get_last_element2)))
  
  constraints <- c()
  for(ii in 1:length(uligands)){
    
    vv1 <- variables$var[which(variables$var_exp==paste0("LR:ligand ", uligands[ii]))]
    vv2 <- variables$var[ind[which(uligands_ctrl==uligands[ii])]]
    
    cc1 <- paste0(vv1, " - ", paste0(vv2, collapse = " - "), " <= 0")
    cc2 <- paste0(vv1, " - ", vv2, " >= 0")
    
    constraints <- c(constraints, unique(c(cc1, cc2)))
    
  }
  
  return(constraints)
  
}