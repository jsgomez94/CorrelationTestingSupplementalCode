###############################################################
###############################################################
###############################################################
##
## Instructions to reproduce testing breast cancer 
##	gene expresión in Section 5.1 of the paper:
## "Effective Permutation Tests for Differences 
##      Across Multiple High-Dimensional Correlation Matrices"
##
## Jose A. Sanchez Gomez, Yongli Zhang, Yufeng Liu
##
###############################################################
###############################################################
###############################################################



###############################################
# 1. Overview of structure
###############################################

Here, we provide a detailed description for the folder structure of our directories, and the functionality of each of the files in the directory.

./breast_gene_cluster/outputs		: Contains the command-line logs for each of the simulation scenarios. Useful for debugging.
./-------------------/req_lib/		: directory containing local package library used for simulations.
./-------------------/result_data/		: directory saving full numerical simulation data in ".csv" format. These results are raw and unprocessed.
./-------------------/src/			: directory containing functions for plotting and data-processing used throughout our analysis.
./-------------------/test_software/	: directory with the file "twosample_test_jose1.R", which defines all the testing functions used throughout our simulation analysis.

./breast_gene_cluster/00_requirements.R		: R code for importing all packages used in our numerical calculations. This code saves the imported packages to the local library "./breast_gene_cluster/req_lib/".
./-------------------/04_Breast_CleanGeneExpression.csv	: csv table containing the processed dataset used for testing.
./-------------------/12_SimulationFunction.R	: R code defining the main numerical experiment function. This function takes an index, imports data, performs tests and saves results in ".csv" format in the "/results_data/" directory.
./-------------------/22_SimulationScript.R	: R script that imports all functions and packages, processes simulation index, runs simulation function and saves results.
./-------------------/30_ClusterPassDebug.sh		: Indexed SLURM cluster calls for numerical experiments for the breast cancer gene expression data. Debugging experiment with 5 simulated replicates.
./-------------------/31_ClusterPassPrelim.sh		: Indexed SLURM cluster calls for numerical experiments for the breast cancer gene expression data. Preliminary simulations with 100 simulated replicates.
./-------------------/32_ClusterPassFull.sh		: Indexed SLURM cluster calls for numerical experiments for the breast cancer gene expression data. Full simulations with 400 simulated replicates.
./-------------------/41_ProcessingResults.R	: R script that takes numerical test results saved in "./breast_gene_cluster/result_data/", averages the results and saves the processed outputs in the files "42_Results_average.csv" for preliminary simulations, and "42_Results_average_400.csv" for full simulations. 
./breast_gene_cluster/43_Plots.R			: R script that uses the file  "42_Results_average_400.csv" to create Figure 4 in our main paper.



###############################################
# 2. Replicating simulations
###############################################

There are two options for generating the plots we provide in our paper: 

(2.A) You can recreate the original plots provided in the paper using the original numerical runs. The data is provided in the files "42_Results_average_400.csv".
        
        Step 1: Copy the provided directory "./breast_gene_cluster/" to your local machine.
        
        Step 2: Open an RStudio session, with "./breast_gene_cluster/" as your working directory. Make sure that the files "42_Results_average.csv" and "42_Results_average_400.csv" are in the directory.

	Step 3: Use this RStudio session to run the file "./breast_gene_cluster/43_Plots.R". This will create the output plot in Figure 5 of the main paper.



(2.B) Rerunning all numerical scenarios using "/32_ClusterPassFull.sh" to regenerate data. Once new data is generated, it must be processed. After processing, you can then, follow (2.A) to generate new simulated plots.

        Step 1: Copy our repository to a LINUX-based computer cluster with SLURM scheduling system. Modify the value of the variable ROOT_PATH in the requirements file "./breast_gene_cluster/req_lib/" to the appropriate directory.

        Step 2: Ensure the following directories exist:
                "./breast_gene_cluster/req_lib/"
                "./breast_gene_cluster/outputs/"
                "./breast_gene_cluster/result_data/"
        
        Step 3: To do an initial debug for verifying that the simulations run properly, open a linux terminal and navigate to "./breast_gene_cluster/". Then, run the line:
                sbatch 30_ClusterPassDebug.sh
		To verify if the process came out without errors, you may check the directory "./breast_gene_cluster/result_data/" to verify if any raw data was created. You can also check the command line logs which might contain error codes in the files "./breast_gene_cluster/outputs/output_***.out". Any bugs or code errors may be reported in these files.

        Step 4: To run preliminary experiments (100 replicates) open a linux terminal and navigate to "./breast_gene_cluster/". Then, run the line:
                sbatch 31_ClusterPassPrelim.sh
		To verify if the process created the simulated data, you may check the directory "./breast_gene_cluster/result_data/". You can also check the command-line logs which might contain error codes in the files "./breast_gene_cluster/outputs/output_***.out". Any bugs or code errors may be reported in these files.

        Step 5: To run full simulations (400 replicates) open a linux terminal and navigate to "./breast_gene_cluster/". Then, run the line:
                sbatch 32_ClusterPassFull.sh
		To verify if the process created the simulated data, you may check the directory "./breast_gene_cluster/result_data/". You can also check the command-line logs which might contain error codes in the files "./breast_gene_cluster/outputs/output_***.out". Any bugs or code errors may be reported in these files.


        Step 6: To obtain plots of the resulting simulations, follow steps outlined in (2.A). 
           