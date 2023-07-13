####################################
## 1. Runnig models script
## Date: 13/087/2023
## Author: Daniel Noble
## Purpose: Run all models in the katana folder
####################################

# Use this script to run all the .R model scripts in the katana folder. This can be done on all models, subsets or individual models depending on the needs of the user. Code would need to be modified to run on a subset or individual models. All models are saved and provided so that this doesn't necessarily need to be run more than once. For full model details, see scripts in the katana folder.

#################
# Full Data Models
#################

# Load the files
	files <- list.files("./katana", pattern = "DA.*\\.R$", full.names = TRUE)[-c(1:2)]

# Run First 8 models as we can have two instanced running in parallel
	for(i in 1:(length(files)/2)){
		source(files[i], verbose = FALSE)
		print(paste0("Finished ", files[i]))
	}

# Second 8 models as a second instance
	for(i in 9:(length(files))){
	source(files[i], , verbose = FALSE)
	print(paste0("Finished ", files[i]))
	}

#################
# Cold Data Models
#################

# Load the files
	files <- list.files("./katana", pattern = "cold.*\\.R$", full.names = TRUE)[-c(1:2)]

# Run First 8 models as we can have two instanced running in parallel
	for(i in 1:(length(files)/2)){
		source(files[i], verbose = FALSE)
		print(paste0("Finished ", files[i]))
	}

# Second 8 models as a second instance
	for(i in 9:(length(files))){
	source(files[i], , verbose = FALSE)
	print(paste0("Finished ", files[i]))
	}

#################
# Hot Data Models
#################

# Load the files
	files <- list.files("./katana", pattern = "hot.*\\.R$", full.names = TRUE)[-c(1:2)]

# Run First 8 models as we can have two instanced running in parallel
	for(i in 1:(length(files)/2)){
		source(files[i], verbose = FALSE)
		print(paste0("Finished ", files[i]))
	}

# Second 8 models as a second instance
	for(i in 9:(length(files))){
	source(files[i], , verbose = FALSE)
	print(paste0("Finished ", files[i]))
	}