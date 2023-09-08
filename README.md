
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