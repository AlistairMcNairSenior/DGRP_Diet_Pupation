
This repository contains data and scripts associated with Havula et al. Genetic variation in macronutrient tolerance in Drosophila melanogaster.

There are two parts to the analysis described here:

Part 1: Estimation of genetic variance.
long.csv contains the data associated with the estimation of genetic variance in response to diet. The analyses are in the script GLMM Tables S2 to S4.

Part 2: GWAS

Scripts are in testing_scripts this folder:

processPhenotypeData.R
- convert Excel spreadsheet into long data format

test_snps_pupation.R
- test SNP associations using two models

test_snps_pupation_downstream.R
- generate figures using these models

These scripts assume the directory structure:

- DGRP_Data
These are genotype tables downloaded from the DGRP version 2

- Data
Phenotype files for pupation and eclosino

- testing_scripts
Scripts for testing associations

- Figures
Output folder for figures

Data are in the Data folder.
