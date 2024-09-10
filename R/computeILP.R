computeILP <- function(variables = variables,
                       background.networks.list = background.networks.list,
                       ccc.input = ccc.input,
                       tf.input = tf.input,
                       ligand.scores = ligand.scores,
                       ccc.prob = ccc.prob,
                       as.input = as.input,
                       solverPath = solverPath,
                       lambda1 = lambda1,
                       lambda2 = lambda2,
                       lambda3 = lambda3,
                       lambda4 = lambda4,
                       mipgap = mipgap,
                       relgap = relgap,
                       populate = populate,
                       nSolutions = nSolutions,
                       timelimit = timelimit,
                       intensity = intensity,
                       replace = replace,
                       threads = threads,
                       condition = condition,
                       save_res = save_res,
                       lambda5 = lambda5){
  
  objective.function <- write_objective_function(variables = variables, 
                                                 background.networks.list = background.networks.list, 
                                                 ccc.input = ccc.input,
                                                 tf.input = tf.input,
                                                 ligand.scores = ligand.scores,
                                                 ccc.prob = ccc.prob,
                                                 lambda1 = lambda1, 
                                                 lambda2 = lambda2, 
                                                 lambda3 = lambda3,
                                                 lambda4 = lambda4,
                                                 lambda5 = lambda5)
  print("Writing of the Objective Function: Done! Now writing constraints...")
  
  print("Writing of the constraints, step 1/11...")
  c1 <- c(write_constraints_1a(variables = variables, 
                               background.networks.list = background.networks.list),
          write_constraints_1b(variables = variables, 
                               background.networks.list = background.networks.list))
  print("Writing of the constraints, step 2/11...")
  c2 <- write_constraints_2(variables = variables, 
                            background.networks.list = background.networks.list)
  print("Writing of the constraints, step 3/11...")
  if(is.null(ccc.input)){
    c3 <- write_constraints_3b(variables = variables, 
                               background.networks.list = background.networks.list)
  } else {
    c3 <- write_constraints_3a(variables = variables, 
                               background.networks.list = background.networks.list)
  }
  print("Writing of the constraints, step 4/11...")
  c4 <- write_constraints_4(variables = variables, 
                            background.networks.list = background.networks.list)
  print("Writing of the constraints, step 5/11...")
  if(is.null(tf.input)){
    c5 <- write_constraints_5b(variables = variables, 
                               background.networks.list = background.networks.list)
  } else {
    c5 <- write_constraints_5a(variables = variables, 
                               background.networks.list = background.networks.list, 
                               tf.input = tf.input)
  }
  print("Writing of the constraints, step 6/11...")
  c6 <- write_constraints_6(variables = variables, 
                            background.networks.list = background.networks.list)
  print("Writing of the constraints, step 7/11...")
  c7 <- write_constraints_7(variables = variables, 
                            background.networks.list = background.networks.list, 
                            as.input = as.input)
  print("Writing of the constraints, step 8/11...")
  c8 <- write_constraints_8(variables = variables, 
                            background.networks.list = background.networks.list)
  print("Writing of the constraints, step 9/11...")
  c9 <- write_constraints_9(variables = variables,
                            background.networks.list = background.networks.list)
  print("Writing of the constraints, step 10/11...")
  c10 <- write_loop_constraints(variables = variables, 
                                background.networks.list = background.networks.list)
  print("Writing of the constraints, step 11/11...")
  c11 <- write_ligand_control_constraints(variables = variables, 
                                          background.networks.list = background.networks.list)
  allC <- unique(c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11))
  # save(allC, file = paste0("all_constraints_", condition, ".RData"))
  
  print("Writing of all the constraints: Done! Now defining the bounds of the variables...")
  bounds <- write_bounds(variables = variables)
  binaries <- write_binaries(variables = variables)
  print("Writing of all the bounds: Done! Now writing the LINDA+ ILP problem...")
  
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
  
  print("Writing of the ILP problem: Done! Now solving the ILP problem with CPLEX.")
  solve_with_cplex(solverPath = solverPath, variables = variables,
                   condition = condition)
  
  print("Problem solving: Done! Now retreiving all the results (final step)...")
  sif <- read_solution_cplex(variables = variables,
                             background.networks.list = background.networks.list,
                             tf.input = tf.input,
                             condition = condition)
  
  cleanupILP(condition = condition)
  print("Results have been retreived and all the auxiliary CPLEX files have been cleared up!")
  
  return(sif)
  
}