computeILP <- function(variables = variables, 
                       background.networks.list = background.networks.list, 
                       tf_scores = tf_scores, 
                       ligand_scores = ligand_scores, 
                       alpha = 10, 
                       beta = 20, 
                       gamma = 0.1,
                       condition = 1,
                       solverPath = solverPath,
                       mipgap = 0.05,
                       relgap = 0.05,
                       intensity = 2,
                       populate = 1000,
                       nSolutions = 100,
                       replace = 1,
                       threads = 0,
                       timelimit = 3600){
  
  objective.function <- write_objective_function(variables = variables, 
                                                 background.networks.list = background.networks.list, 
                                                 tf_scores = tf_scores, 
                                                 ligand_scores = ligand_scores, 
                                                 alpha = alpha, 
                                                 beta = beta, 
                                                 gamma = gamma)
  
  c1 <- write_constraints_1(variables = variables, 
                            background.networks.list = background.networks.list)
  c2 <- write_constraints_2(variables = variables, 
                            background.networks.list = background.networks.list)
  c3 <- write_constraints_3(variables = variables, 
                            background.networks.list = background.networks.list)
  c4 <- write_constraints_4(variables = variables, 
                            background.networks.list = background.networks.list)
  c5 <- write_constraints_5(variables = variables, 
                            background.networks.list = background.networks.list, 
                            tf_scores = tf_scores)
  c6 <- write_loop_constraints(variables = variables, 
                               background.networks.list = background.networks.list)
  allC <- unique(c(c1, c2, c3, c4, c5, c6))
  
  bounds <- write_bounds(variables = variables)
  binaries <- write_binaries(variables = variables)
  
  # write the .lp file
  data = paste0("testFile_", condition, ".lp")
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  write(objective.function, data, append = TRUE)
  write("Subject To", data, append = TRUE)
  write(allC, data, append = TRUE)
  write("Bounds", data, append = TRUE)
  write(bounds, data, append = TRUE)
  write("Binaries", data, append = TRUE)
  write(binaries, data, append = TRUE)
  write("End", data, append = TRUE)
  
  # write cplexCommand file
  data2 = paste0("cplexCommand_", condition, ".txt")
  write(paste0("read testFile_", condition, ".lp"), data2)
  write(paste0("set mip tolerances mipgap ", mipgap), data2, append = TRUE)
  write(paste0("set mip pool relgap ", relgap), data2, append = TRUE)
  write(paste0("set mip pool replace ", replace), data2, append = TRUE)
  write(paste0("set mip limits populate ", populate), data2, append = TRUE)
  write(paste0("set mip pool capacity ", nSolutions), data2, append = TRUE)
  write(paste0("set mip pool intensity ", intensity), data2, append = TRUE)
  write(paste0("set timelimit ", timelimit), data2, append = TRUE)
  write(paste0("set threads ", threads), data2, append = TRUE)
  write("populate", data2, append = TRUE)
  write(paste0("write results_", condition, ".txt sol all"), data2,
        append = TRUE)
  write("quit", data2, append = TRUE)
  
  solve_with_cplex(solverPath = solverPath, variables = variables,
                   condition = condition)
  
  sif <- read_solution_cplex(variables = variables,
                             background.network = background.networks.list,
                             condition = condition)
  
  cleanupILP(condition = condition)
  
  return(sif)
  
}