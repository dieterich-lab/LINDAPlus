library(LINDA)
library(XML)

# get input data
load(file = system.file("extdata", "toy.background.networks.list.RData", package = "LINDAPlus"))
load(file = system.file("extdata", "toy.ligand.scores.RData", package = "LINDAPlus"))
load(file = system.file("extdata", "toy.tf.scores.RData", package = "LINDAPlus"))
load(file = system.file("extdata", "toy.top.tf.RData", package = "LINDAPlus"))

# get expected result
load(file = system.file("result_expected.RData",
                        package="LINDA"))

# obtain actual reesult
result_actual = runLINDAPlus(background.networks.list = background.networks.list, 
                             ligand.scores = ligand_scores, 
                             tf.scores = tf_scores, 
                             solverPath = "~/Downloads/cplex", 
                             top = top)

#testing
test_that("Comparison of the results", {
  expect_equal(result_actual, result_expected)
})