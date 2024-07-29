#' Performs the LINDA analysis once the inputs are given.
#'
#'@param background.networks.list The list consists of two elements:
#'a)background.networks: This should be a named (cell-types) list containing
#'data-frames joint PPI-DDI's for each cell-type. Each data-frame represents
#'the set of interaction knowledge for each cell-type and it should contain at
#'least 4 columns with the following ID's: 'pfam_source', 'pfam_target',
#'gene_source' and 'gene_target'. In the case where we have no name for a
#'specific domain or where this is not applicable, please set the corresponding
#'values in the 'pfam_source' or 'pfam_target' as NA's.
#'b)ligand.receptors: This also should be a named list ('ligands' and Receptors') 
#'which contains character vectors where the set of elements that 
#'corresponds to Ligands and Receptors have been defined as such.
#'
#'
#'@param ligand.scores Users can provide information about the abundance of 
#'ligands in the extra-cellular space as made evident by Secretomics data 
#'through a data-frame object. More abundant ligands/extra-cellular molecules 
#'are more likely to initiate conformational changes in receptors. The 
#'data-frame provided should contain two columns: 'ligands' (providing the 
#'ligand ID's) and 'score' (providing the score associated to each ligand, i.e. 
#'abundance). The higher the score of the ligand, the more likely it will be for
#' a ligand to appear in the solution.
#'
#'
#'@param tf.scores This information should be provided as a named list (for each 
#'cell-type) and which contains data-frames indicating the enrichment scores for 
#'each TF at each cell-type. The data-frames should contain at least two 
#'columns: 'tf' (indicating the TF ID) and the numerical 'score' (indicating the
#' enrichment scores for each TF).
#' 
#' 
#' #'@param top (optional) Users can provide a named (also by cell-type) 
#' numerical vector to indicate the number of TF's to consider as significantly 
#' regulated based on their absolute enrichment values. In case that this 
#' parameter has not been defined, then by default all the TF's provided in the 
#' data-frames list will be considered as significantly regulated.
#'
#'
#'@param solverPath (optional) location path to the desired solver. By default,
#'the path to the cplex solver is set to: solverPath = "/usr/bin/cplex".
#'
#'
#'@param alpha the penalization term of the primary objective of the objective
#'function - TF inclusion. This penalty factor is suggested to be set to a
#'higher value compared to other penalty parameters in order to strongly
#'penalize the inclusion of not signifcantly regulated TF's. By default,
#'alpha=10.
#'
#'
#'@param beta the penalization terms of the secondary objective of the
#'objective function - size penalty. The aim of this objective term is to
#'penalize the inclusion of spurious DDI's in the final solution. By default,
#'beta=5.
#'
#'
#'@param gamma the penalization terms of the secondary objective of the
#'objective function - size penalty. The aim of this objective term is to
#'penalize the inclusion of spurious DDI's in the final solution. By default,
#'gamma=0.1.
#'
#'
#'@param mipgap CPLEX parameter which sets an absolute tolerance on the gap
#'between the best integer objective and the objective of the best node
#'remaining. When this difference falls below the value of this parameter, the
#'mixed integer optimization is stopped. By default, mipgap=0.
#'
#'
#'@param relgap CPLEX parameter which sets a relative tolerance on the objective
#'value for the solutions in the solution pool. Solutions that are worse (either
#'greater in the case of a minimization, or less in the case of a maximization)
#'than the incumbent solution by this measure are not kept in the solution pool.
#'For example, if relgap=0.001 (or 0.1 percent), meaning that then solutions
#'worse than the incumbent by 0.1 percent or more will be discarded. By default,
#'relgap=0.
#'
#'
#'@param populate CPLEX parameter which sets the maximum number of mixed integer
#'programming (MIP) solutions generated for the solution pool during each call
#'to the populate procedure. Populate stops when it has generated the amount of
#'solutions set in this parameter. By default, populate=500.
#'
#'
#'@param nSolutions the number of solutions to be provided by LINDA. By default,
#'nSolutions=100.
#'
#'
#'@param timelimit CPLEX parameter which sets the maximum optimization time in
#'seconds. By default, timelimit=3600.
#'
#'
#'@param intensity CPLEX parameter which controls the trade-off between the
#'number of solutions generated for the solution pool and the amount of time or
#'memory consumed. Values from 1 to 4 invoke increasing effort to find larger
#'numbers of solutions. Higher values are more expensive in terms of time and
#'memory but are likely to yield more solutions. By default, intensity=0 (let
#'CPLEX choose). By default, intensity=1.
#'
#'
#'@param replace CPLEX parameter which designates the strategy for replacing a
#'solution in the solution pool when the solution pool has reached its capacity.
#'The value 0 replaces solutions according to a first-in, first-out policy. The
#'value 1 keeps the solutions with the best objective values. The value 2
#'replaces solutions in order to build a set of diverse solutions. By default,
#'replace=1.
#'
#'
#'@param threads CPLEX parameter which manage the number of threads that CPLEX
#'uses. By default, threads=0 (let CPLEX decide). The number of threads that
#'CPLEX uses is no more than the number of CPU cores available on the computer
#'where CPLEX is running.
#'
#'
#'@param condition a parameter which can be used in the case when LINDA is
#'desirde to run over multiple analyses in parallel. It is useful to distinguish
#'between the multiple ILP problems defined for each case as well as the
#'solutions obtained. By default, conditions=1.
#'
#'
#' @return Results list containing each of the LINDA unique solutions as well as
#' the combined ones.
#'
#' @examples
#' library(LINDAPlus)
#' library(XML)
#'
#' load(file = system.file("extdata", "toy.background.networks.list.RData", package = "LINDAPlus"))
#' load(file = system.file("extdata", "toy.tf.scores.RData", package = "LINDAPlus"))
#' load(file = system.file("extdata", "toy.top.tf.RData", package = "LINDAPlus"))
#'
#' res <- runLINDA(background.networks.list = background.networks.list,
#' ligand.scores = ligand_scores,
#' tf.scores = tf_scores,
#' solverPath = "~/Downloads/cplex",
#' top = top)
#'
#'
#' @export

runLINDAPlus <- function(background.networks.list = background.networks.list,
                         tf.scores = tf.scores,
                         lr.scores = NULL,
                         ligand.scores = NULL,
                         ccc.scores = NULL,
                         as.input = NULL,
                         solverPath = "/usr/bin/cplex",
                         top.tf = NULL,
                         lambda1 = 10,
                         lambda2 = 15,
                         lambda3 = 1,
                         lambda4 = 0.1,
                         mipgap = 0,
                         relgap = 0,
                         populate = 500,
                         nSolutions = 100,
                         timelimit = 3600,
                         intensity = 1,
                         replace = 1,
                         threads = 0,
                         condition = 1,
                         save_res = FALSE){
  
  options(scipen=999)
  
  print("LINDA analysis start !!")
  print("Now checking all the inputs")
  
  all_inputs <- checkInputs(background.networks.list = background.networks.list,
                            ligand.scores = ligand.scores,
                            tf.scores = tf.scores,
                            lr.scores = lr.scores,
                            solverPath = solverPath,
                            top.tf = top.tf, 
                            ccc.scores = ccc.scores,
                            as.input = as.input,
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
                            save_res = save_res)
  
  
  background.networks.list <- all_inputs$background.networks.list
  tf.scores <- all_inputs$tf.scores
  lr.scores <- all_inputs$lr.scores
  ligand.scores <- all_inputs$ligand.scores
  ccc.scores <- all_inputs$ccc.scores
  as.input <- all_inputs$as.input
  solverPath <- all_inputs$solverPath
  top.tf <- all_inputs$top.tf
  lambda1 <- all_inputs$lambda1
  lambda2 <- all_inputs$lambda2
  lambda3 <- all_inputs$lambda3
  lambda4 <- all_inputs$lambda4
  mipgap <- all_inputs$mipgap
  relgap <- all_inputs$relgap
  populate <- all_inputs$populate
  nSolutions <- all_inputs$nSolutions
  timelimit <- all_inputs$timelimit
  intensity <- all_inputs$intensity
  replace <- all_inputs$replace
  threads <- all_inputs$threads
  condition <- all_inputs$condition
  save_res <- all_inputs$save_res
  
  print("Checking of all the inputs: Done! Now processing background network...")
  
  
  background.networks.list <- process_background_network(background.networks.list = background.networks.list, 
                                                         tf.scores = tf.scores, lr.scores = lr.scores, ccc.scores = ccc.scores)
  
  print("Processing the background network: Done! Now creating all the variables...")
  
  variables <- create_variables(background.networks.list = background.networks.list)
  
  print("Writing of all the variables: Done! Now writing the objective function...")
  
  res <- computeILP(variables = variables, background.networks.list = background.networks.list, 
                    tf.scores = tf.scores, ligand.scores = ligand.scores, lr.scores = lr.scores, 
                    ccc.scores = ccc.scores, as.input = as.input, lambda1 = lambda1, 
                    lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4, condition = condition, 
                    solverPath = solverPath, mipgap = mipgap, relgap = relgap, intensity = intensity, 
                    populate = populate, nSolutions = nSolutions, replace = replace, threads = threads, 
                    timelimit = timelimit)
  
  if(save_res){
    save(res, file = paste0("res_", condition, ".RData"))
  }
  
  return(res)
  
}