#################################################
##
## In the following document, we introduce the
## script that has to be run to get the 
## simulation outputs.
##
#################################################
#################################################
## Step 1: import tools.
source("00_requirements.R")
source("test_software/twosample_test_jose1.R")
source("12_SimulationFunction.R")

#################################################
#################################################
## Step 2: proces inputs.
input <- commandArgs(trailingOnly = TRUE)
run_id <- as.numeric(input[1])

#################################################
#################################################
## Step 3: Simulations
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#------------------Running simulations")

full_simulation(run_id = run_id)
