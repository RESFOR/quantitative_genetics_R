# quantitative_genetics_R

# multi-trait and milti-site ABLUP and HBLUP Copyright @ Dr. Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar).

# The multi-trait and milti-site R-script fits the standard pedigree-based ABLUP and combined pedigree-genomic based HBLUP models using the breedR R-package. The outputs of this R-script are the genetic and residual (co)variance components, and the breeding values and its accuracies for both models. 

# As an example, the RES-FOR lodgepole pine dataset (data_multi-trait_multi-site.txt) with 3 traits (DBH30, HT30, NSWGR36) and 4 sites (JUDY, VIRG, SWAN, and TIME) and the H-matrix calculated in blupf90 programs using the 8K ancestry selected SNPs.

# data_multi-trait_multi-site.txt is the demo dataset for the multi-trait and multi-site R script.

# White_Spruce_Phenotype_Pedigree_PLoS2022.txt is the tree pedigree and phenotype data published in Cappa et al. PLoS 2020.

# Spatial-competition-correction-mixed-model.R Copyright @ Dr. Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar).

# This R-script fit the standard, and the spatial (spline and autoregressive -AR1-) and/or competition individual-tree mixed models using the breedR R-package. Moreover, this R-script presents several modalities for diagnosing spatial and competition problems in forest genetic trials. Finally, the last lines of this R-script calculats the Design and Spatial phenotype adjustments. We also provide the total height at age 30 from the RESFOR spruce Red Earth (REDE) site as an example (data_spatial_competition_correction.txt).

# data_spatial_competition_correction.txt is the demo dataset for the Spatial-competition-correction-mixed-model.R.

# AIM.coefficient.R performs SNP variable selection based on the AIM (ancestry informative marker, Rosenberg et al. 2003).

# AIM.coefficient.R copyright @ Blaise Ratcliffe (b.ratcliffe@gmail.com).

# Multiple-site Single-trait GBLUP.R performs mixed models using breedR R package. Copyright @ Dr. Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar).

# The G_matrix_RESFOR_lodgepole_pine.RDS is the genomic relationship matrix of 1,490 lodgepole pine trees genotyped by 25,099 GBS SNPs. Copyright @ Dr. Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar). The G matrix file is too large to be visible on GitHub. The file is downloadable.  

# The G matrix was used in the R script that fits a genomic-based GBLUP model using the breedR R-package. The output of this R-script includes the genetic and residual (co)variance components (and function of them, heritability estimates and genetic correlations between sites -GxE-), and the breeding values and its´ standard errors from the GBLUP model.

# R-script_EDA.R is a R-scrip for Exploratory Data Analysis (EDA). The EDA is an important step before start any model fit, which may help to identify errors, detect outliers and any anomalous events in the data, and find relations between the phenotypic variables. This R-script include univariate and multivariate visualization and summary statistics of continuous and categorical phenotypic traits. We have also included a preliminary analysis of spatial diagnosis. In this case, as example, I used the RAW RESFOR lodgepole pine dataset with three traits (DBH30, HT30, and WGR36) and the 4 sites (JUDY, VIRG, SWAN, and TIME). The data (data_EDA_analysis.txt) used here is also available. 

# R-script_EDA.R Copyright @ Dr. Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar).

# data_EDA_analysis.txt is the example data set used in the R-script_EDA.R. Copyright @ Dr. Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar).
