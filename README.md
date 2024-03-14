# LINDAPlus
LINDA+ R-package - Enio GJERGA

This repository contains the scripts for the ILP (Integer Linear Programming) implementation of the LINDA+ R package and accompanying scripts that implement the method. ILP is a mathematical optimisation algorithm in which the objective function and constraints are linear and the variables are integers.

### License

Distributed under the GNU GPLv3 License.

### Installation

#### 1. Solver Prerequisites
Before installing LINDA+, please keep in mind the following solver pre-requisites:

LINDA+ requires the interactive version of IBM Cplex as an ILP problem optimiser. The IBM ILOG Cplex is freely available through [Academic Initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB).

Details about how to obtain the free full license under the Academic Initiative have been provided [here](https://community.ibm.com/community/user/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students).

Once the CPLEX solver has been acquired by the user, they must save the executable files in any desired location in the machine they are using and then they can run. Then, thy should make the acquired _cplex_ solver executable, but running the following (in MAC/Unix systems): ```chmod +x cplex```. Finally, you can perform the LINDA+ analyses after specifying the path to the executable file (through the *solverPath* parameter, see examples below).

#### 2. Package Depedencies
Additionally before installation, the users must install the following R-package depedencies:
[igraph](https://igraph.org/r/), [aggregation](https://cran.r-project.org/web/packages/aggregation/index.html) and
[XML](https://cran.r-project.org/web/packages/XML/index.html).

#### 3. Package installation
Once the required solvers have been obtained and the mentioned R-package depedencies have been installed, then the users can proceed with the installation of LINDA.

Currently users can install LINDA directly from the source after downloading the source (tar file) and typing in ```R``` command line the following:

```R
# directly from GitHub
devtools::install_github("dieterich-lab/LINDAPlus", build_vignettes = FALSE)
```

```R
# or download the source file from GitHub and install from source
devtools::install_local(path = "path_to_extracted_LINDAPlus_directory", build_manual = TRUE, build_vignettes = TRUE, force = TRUE)
```

**NOTE:** If you wish for the the Vignettes to be built and for the test example to run successfully, please put the _cplex_ executable file under the /Downloads/ directory and **only** then you can set ```build_vignettes = TRUE```. In the case when building the vignettes is not possible, you can access it [here](https://github.com/dieterich-lab/LINDAPlus/blob/main/vignettes/LINDAPlus.html).

Upon installation of the package, you can check on the package Vignettes:
```R
# check vignettes
vignette("LINDAPlus")

# check documentation of the main function
?runLINDAPlus
```

## Running LINDA+

The LINDA library can be initialized by:

```R
library(LINDAPlus)
```

A documentation of the current version of the main _runLINDAPlus()_ function can be obtained by simply typing ```?runLINDAPlus``` in R once the package has been installed.

## Examples
Below we provide the codes about running a test toy example (assuming that the _cplex_ solver has been stored in ```"~/Downloads/cplex"```) after calling for the _lindaPlus_ package by using the default optimization parameter settings of the _runLINDAPlus()_ function:

```R
load(file = system.file("extdata", "toy.background.networks.list.RData", package = "LINDAPlus")) #backgroud network for each cell-type + ligands&receptors objects 
load(file = system.file("extdata", "toy.ligand.scores.RData", package = "LINDAPlus")) #ligand scores
load(file = system.file("extdata", "toy.tf.scores.RData", package = "LINDAPlus")) #TF scores for each cell-type
load(file = system.file("extdata", "toy.top.tf.RData", package = "LINDAPlus")) #The number of TF's to consider as top-regulated for each cell-type

res <- runLINDAPlus(background.networks.list = background.networks.list, 
                    ligand.scores = ligand_scores, 
                    tf.scores = tf_scores, 
                    solverPath = "~/Downloads/cplex", 
                    top = top)

# Show the combined table of inferred inter-cellular Ligand-Receptors interactions
res$combined_solutions$ligand_receptors

# Show the combined table of inferred intra-cellular PPI's for CellA 
res$combined_solutions$CellA

# Show the combined table of inferred intra-cellular PPI's for CellB
res$combined_solutions$CellB
```

Alternatively users can expand on the gap tolerance and/or the diversity of solutions to report in order to obtain a wider spectrum of alternative solutions.
```R
load(file = system.file("extdata", "toy.background.networks.list.RData", package = "LINDAPlus")) #backgroud network for each cell-type + ligands&receptors objects 
load(file = system.file("extdata", "toy.ligand.scores.RData", package = "LINDAPlus")) #ligand scores
load(file = system.file("extdata", "toy.tf.scores.RData", package = "LINDAPlus")) #TF scores for each cell-type
load(file = system.file("extdata", "toy.top.tf.RData", package = "LINDAPlus")) #The number of TF's to consider as top-regulated for each cell-type

res <- runLINDAPlus(background.networks.list = background.networks.list, 
                     ligand.scores = ligand_scores, 
                     tf.scores = tf_scores, 
                     solverPath = "~/Downloads/cplex", 
                     top = top, 
                     mipgap = 0.01, 
                     relgap = 0.01, 
                     replace = 2)

# Show the combined table of inferred inter-cellular Ligand-Receptors interactions
res$combined_solutions$ligand_receptors

# Show the combined table of inferred intra-cellular PPI's for CellA 
res$combined_solutions$CellA

# Show the combined table of inferred intra-cellular PPI's for CellB
res$combined_solutions$CellB
```
