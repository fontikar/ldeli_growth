
# Heritability and developmental plasticity of growth in an oviparous lizard

This repository contains the data, code and models associated with the following manuscript:

Fonti Kar, Shinichi Nakagawa, Daniel W.A. Noble. Heritability and developmental plasticity of growth in an oviparous lizard. EcoEvoRxiv. https://doi.org/10.32942/X2PP47. 

## Relevant code files

The key scripts needed to reproduce the results are the following:
1) `Results_brms.qmd`: This script contains the code to run the models and generate the figures and tables in the main manuscript.
2) `ESM_brms.qmd`: This script contains the code to run the models and generate the figures and tables in the electronic supplementary materials.
3) `R/Analsysis/functions.R`: This file contains all the functions used for processing models and making calculations
4) `R/run_mods.R`: Should run all the model files in the `Katana` folder. The models are then saved to the `output/rds` folder. 

There are also a number of other scripts that are used clean the data and generate G matrices which are saved in the `Analysis/ Data clean and process` folder. The various files above will load the relevant data and model files. 

## Models

All models are run in a Bayesian framework. As such, they take *a lot* of time to fit. We do not recommend rerunning the scripts in the `Katana` folder. Instead, we recommend using the `output/rds` files which contain the fitted models. These are loaded already within the code files mentioned above.

If you need to know more about what each model is doing, please see the `models_table.xlsx` file which describes the model name and structure of the model. 

## SNP Data

Unfortauntely, the SNP genotyping data is astronomically large! Which means that these data cannot be stored on GitHub (sorry!). If you would like access to these data, please contact the corresponding author (Daniel Noble) and we can arrange to send you the data or see the [OSF site](https://osf.io/hjkxd/?view_only=12a6b6010567474fac9fecd54472aa3d). This contains a folder called `dartR` which contains the needed files for loading G matrices and SNP data. Please place this folder in the `output` folder of the GitHub repository and the path names should match if your working directory is set to the root of the GitHub repository.