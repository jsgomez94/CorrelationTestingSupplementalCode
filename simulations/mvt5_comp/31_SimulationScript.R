#################################################
##
## In the following document, we introduce the
## script that has to be run to get the 
## simulation outputs.
##
#################################################

source("1_requirements.R")

source("mvt5_comp/11_twosample_test_jose1.R")
source("mvt5_comp/20_model5.R")
source("mvt5_comp/22_CreatingParameters.R")
source("mvt5_comp/23_SimulationFunction.R")

#################################################
#################################################
## Step 1: Read imputs:
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#-----------------------Reading inputs")

input <- commandArgs(trailingOnly = TRUE)
run_id <- as.numeric(input[1])
run_type <- as.numeric(input[2])
id_task <- (run_id %/% 100)
run_id  <-  (run_id %% 100)

#plist <- list(
#    p15 = 0:71, 
#    p30 = 72:143,
#    p60 = 144:215,
#    p120 = 216:287,
#    p240 = 288:359,
#    p480 = 360:431,
#    p960 = 432:503)

saving_folder <- ifelse(run_type == 1, "mvt5_comp/txt_exps/", "mvt5_comp/txt_full/")
output_folder <- ifelse(run_type == 1, "mvt5_comp/output_exps/", "mvt5_comp/output_full/")
if (!dir.exists(saving_folder)) {
       dir.create(saving_folder)
}
if (!dir.exists(output_folder)) {
       dir.create(output_folder)
}

#################################################
#################################################
## Step 2: Run simulation scenarios:
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#------------------Running simulations")
assign(paste0("output", run_id),
       full_simulation(id_task, run_id, run_type))
print(get(paste0("output", run_id)))