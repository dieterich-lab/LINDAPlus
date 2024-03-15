## Overview

[LINDA+](https://github.com/dieterich-lab/LINDAPlus), is an R-Package used to identify functional intra-cellular Protein-Protein (PPI) and Domain-Domain (DDI) Interactions as well as inter-cellular Ligand-Recpetor (LR) Interactions between various types of cells based on sc-RNAseq and Secretomics data. It is an upgrade to the previous [LINDA](https://dieterich-lab.github.io/LINDA/) (Linear Integer programming for Network reconstruction using transcriptomics and Differential splicing data Analysis) methodology which inferred mechanistic upstream regulatory signalling processes that drive gene expression changes by also taking into account the contribution of alternative splicing events to signal protein variability and information flow modulation within a single cell from bulk data.

### Installation and Requirements

LINDA+ can be installed either locally by downloading its [source code](https://github.com/dieterich-lab/LINDAPlus) or directly from GitHub via [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html).


#### Install through download
1.  Download LINDA package (*Main* or *Development* branch) by clicking on **Code** and then **Download ZIP**.
2.  Start R.
3.  Unzip package and set working directory to where LINDA+ has been downloaded in the R workspace.
4.  Install LINDA by running: `install.packages('LINDAPlus/', repos = NULL, type="source")`.
5. Alternatively you can install locally with [devtools](https://cran.r-project.org/web/packages/devtools/index.html) by running: `devtools::install_local(path = "path_to_extracted_LINDAPlus_directory", build_manual = TRUE, build_vignettes = FALSE, force = TRUE)`.

**NOTE:** If you wish for the the Vignettes to be built and for the test example to run successfully, please put the cplex executable file under the _/Downloads/_ directory and only then you can set `build_vignettes = TRUE`. In the case when building the vignettes is not possible, you can access it within the LINDA+ [R-package](https://github.com/dieterich-lab/LINDAPlus) [here](https://github.com/dieterich-lab/LINDAPlus/blob/master/vignettes/LINDAPlus.html).

#### Install from GitHub
1.  Make sure [devtools](https://cran.r-project.org/web/packages/devtools/index.html) is installed.
2.  Start R.
3.  Run `install_github("dieterich-lab/LINDAPlus", build_vignettes = TRUE)`. (**NOTE:** For the test example to run successfully, please put the _cplex_ executable file under the _/Downloads/_ directory. Otherwise set `build_vignettes = FALSE`).

#### Prerequistes
Upon installation, it is required to obtain the interactive version of IBM Cplex as an ILP problem optimiser. The IBM ILOG Cplex is freely available through [Academic Initiative](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_%7C470%7C135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB). Details about how to obtain the free full license under the Academic Initiative have been provided [here](https://community.ibm.com/community/user/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students).

### Documentation
The LINDA pipeline and it prerequisites have been documented on its Vignettes. Also documentation about to test cases (one Toy example and one Real-case application) have been provided.
1.  For an explanation of the LINDA pipeline and a to run a small Toy test case-study, you can check the vignettes by simply running: `vignette("LINDAPlus")`.
2.  Also a step-by-step explanation of LINDA application over a real case-study has been provided in: To be updated ...

### Support or Contact

Having trouble with LINDA? You can open an Issue on the dedicated [LINDA+ repository](https://github.com/dieterich-lab/LINDAPlus), or simply drop us an email on **E.Gjerga@uni-heidelberg.de** & **Christoph.Dieterich@uni-heidelberg.de** and we’ll help you sort it out.
