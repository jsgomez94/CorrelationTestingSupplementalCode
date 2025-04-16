# HubEstimationCodeSubmission

This folder contains the code corresponding to the numerical simulations results found in Section 5 in the main paper, and Sections S2-S8 of our Supplementary Materials. First, we here provide a description of the overall structure of the files. Second, we describe in detail the steps required to replicate the simulations, 

## Overall Structure of Files/Directories:

- Directory <code>./req_lib</code>: directory where the required packages will be saved as local package installations. 

- Directory <code>./mvn_comp</code>: directory that contains all the files for reproducing the comparative simulations for correlation/covariance testing with Gaussian distributions found in Section 4 of the main paper, and Sections S2, S3, S4 of our Supplementary Materials. 

- Directory <code>./mvn_supps</code>: directory that contains all the files for reproducing the simulations exploring the performance of different choices of entrywise transformations, matrix norm, and number of permutation resampling B found in Sections S5, S6, and S7 of our Supplementary Materials. 

- Directory <code>./mvt5_comp</code>: directory that contains all the files for reproducing the comparative simulations for correlation/covariance testing with multivariate T5 distributions found in Section S8 of our Supplementary Materials. 

- File <code>./1_requirements.R</code>: file that creates local installations of all the packages necessary for the simulation experiments to the directory <code>./req_lib</code>.

- Files numbered <code>./2-5***.R</code>: files that create the numerical result graphical summaries for each of our simulation scenarios. Each file is named according to which simulation scenario setting it processes and visualizes. 


## Simulation Rerunning Instructions:

The simulations are run in a linux cluster with SLURM scheduling system. 


0. **Clone Repository:** First, clone the github repository <em>CorrTestingCodeSubmission</em> to a linux cluster with SLURM scheduling system.

1. **Install Requirements:** Install all packages required for our simulations by running

    - <code>Rscript 1_requirements.R</code>

    Ensure that all packages are properly installed by reviewing the directory <code>./req_lib/</code>.

2. **Debugging to run mvn_comp:** Due to the large volume of simulations, and potential emerging errors when cloning the repository into a new environment, it may be necessary to perform some debugging rounds and verify that the code runs properly. 

    In the linux cluster terminal, enter the directory <code>./simulations/</code>. To ensure that things run smoothly, run the command,

    - <code>sbatch mvn_comp/40_ClusterPassExps_debug.sh</code>

    This generates a small simulation run, which will help verify if the simulation code is successfully running and saving outputs. To ensure the code ended successfully, verify the directory <code>./mvn_comp/txt_exps/</code>, which should contain a total of 200 files with simulation outputs in TXT format. In case of error, such as no TXT files in <code>./mvn_comp/txt_exps/</code> or corrupted files, explore the directory <code>./mvn_comp/output_exps/</code>, which should contain 100 files with name format <code>output_p15_**.out</code> numbered from 0-99, and contain the log file outputs of each of the runs. Use these command-line log files to help your debugging efforts.

    Just as with the directory <code>./mvn_comp/</code>, similar debugging steps can be taken to verify whether simulations can be run for <code>./mvn_supps/</code> and <code>./mvt5_comp/</code> by running the following lines of code:

    - <code>sbatch mvn_supps/40_ClusterPassExps_debug.sh</code>
    - <code>sbatch mvt5_comp/40_ClusterPassExps_debug.sh</code>

    Then, exploring the TXT numerical data files found in <code>./mvn_supps/txt_exps/</code> and <code>./mvt5_comp/txt_exps/</code> can help verify that the code finished successfully. If files are missing or are corrupted, we recommend verifying the log files in <code>./mvn_supps/output_exps/</code> and <code>./mvt5_comp/output_exps/</code> to determine possible causes of error.


3. **Running Simulation Experiments:** Once debugging experiments ran succesfully, you can run preliminary simulation experiments. We have a total of 168 simulation scenarios, corresponding to different choices of dimension p, sample size n, underlying hypotheses (0 for null,1 for alternative sparse, 2 for alternative dense), and correlation matrix model (1,2,...,5). For each of these simulation scenarios, we perform 200 simulation replicates. To do this, we request 100 cluster nodes, and ask each to perform 2 simulation replicates per simulation scenarios. To run this, in a cluster terminal, enter the directory <code>./simulations/</code>, and run the lines,

    - <code>sbatch mvn_comp/41_ClusterPassExps_p15.sh</code>
    - <code>sbatch mvn_comp/42_ClusterPassExps_p30.sh</code>
    - <code>sbatch mvn_comp/43_ClusterPassExps_p60.sh</code>
    - <code>sbatch mvn_comp/44_ClusterPassExps_p120.sh</code>
    - <code>sbatch mvn_comp/45_ClusterPassExps_p240.sh</code>
    - <code>sbatch mvn_comp/46_ClusterPassExps_p480.sh</code>
    - <code>sbatch mvn_comp/47_ClusterPassExps_p960.sh</code>
    
    The logs for the Correlation Testing reduced experiments are saved in <code>./mvn_comp/output_exps/</code>, indexed from 100 to 16899. The simulation data is saved in TXT files in <code>./mvn_comp/txt_exps/</code>. 
    
    Similar log and TXT can be found for the simulations of the directories <code>./mvn_supps/</code>  and <code>./mvt5_comp/</code>. To run preliminary simulation experiments, simply run the commands:

    - <code>sbatch mvn_supps/41_ClusterPassExps_p15.sh</code>
    - <code>sbatch mvn_supps/42_ClusterPassExps_p30.sh</code>
    - <code>sbatch mvn_supps/43_ClusterPassExps_p60.sh</code>
    - and so on...

    - <code>sbatch mvt5_comp/41_ClusterPassExps_p15.sh</code>
    - <code>sbatch mvt5_comp/42_ClusterPassExps_p30.sh</code>
    - <code>sbatch mvt5_comp/43_ClusterPassExps_p60.sh</code>
    - and so on...

    

4. **Generating Simulation Experiment Plots:** Once all preliminary simulation experiments are complete and the data is saved, you can generate preliminary simulation plots. For this, in the command line terminal enter the directory <code>./simulations/</code>, and run:

    - <code>Rscript 2_twosample_plots_sizepower_mvn_comparison.r</code> 
    
    This saves aggregated data in the directory <code>./results/mvn_comp_exp/</code>. Plots that verify performance in terms of empirical rejection rate are then saved in the directory <code>./results/mvn_comp/Method_comparison/</code>. 
    
    To generate the aggregated data and figures corresponding to the supplementary comparisons of Sections S5, S6, and S7 of our Supplementary Materials, simply enter the directory <code>./simulations/</code> in the command line terminal and run:

    - <code>Rscript 3_twosample_plots_sizepower_mvn_supps_Bcomp.r</code> 
    - <code>Rscript 4_twosample_plots_sizepower_mvn_supps_Norm.r</code> 
    
    Finally, to generate the aggregated data and figures corresponding to the supplementary comparisons of Section S8 of our Supplementary Materials, simply enter the directory <code>./simulations/</code> in the command line terminal and run:

    - <code>Rscript 5_twosample_plots_sizepower_mvt5_comparison.r</code> 

    **NOTE 1:** In all the figure-generating scripts described above, there is a variable <code>sim_type</code> which needs to be adjusted depending on whether the simulation is a reduced experiment or a full simulation run. Ensure you modify this value accordingly to obtain the appropriate results. 

    **NOTE 2:** The object <code>sim_type</code> is defined more than once within each of the plotting files. Before you run, ensure that you have set the value of <code>sim_type</code> to <code>sim_type <- "exps"</code> in all instances before running. 



5. **Running Full Simulations:** You can also run full numerical simulations. For each of the 168 simulation scenarios, we perform 2000 simulation replicates. To do this, we request 100 cluster nodes, and ask each to perform 20 simulation replicates per simulation scenario. To run this, in a cluster terminal, enter the directory <code>./simulations_4/</code>, and run the lines,
    
    - <code>sbatch mvn_comp/51_ClusterPassFull_p15.sh</code>
    - <code>sbatch mvn_comp/52_ClusterPassFull_p30.sh</code>
    - <code>sbatch mvn_comp/53_ClusterPassFull_p60.sh</code>
    - <code>sbatch mvn_comp/54_ClusterPassFull_p120.sh</code>
    - <code>sbatch mvn_comp/55_ClusterPassFull_p240.sh</code>
    - <code>sbatch mvn_comp/56_ClusterPassFull_p480.sh</code>
    - <code>sbatch mvn_comp/57_ClusterPassFull_p960.sh</code>
    
    The logs for the Correlation Testing full simulation runs are saved in <code>./mvn_comp/output_full/</code>, indexed from 100 to 16899. The simulation data is saved in TXT files in <code>./mvn_comp/txt_full/</code>. 

    Similar log and TXT can be found for the simulations of the directories <code>./mvn_supps/</code>  and <code>./mvt5_comp/</code>. To run full simulation runs with 2000 replicates, simply run the commands:

    - <code>sbatch mvn_supps/51_ClusterPassFull_p15.sh</code>
    - <code>sbatch mvn_supps/52_ClusterPassFull_p30.sh</code>
    - <code>sbatch mvn_supps/53_ClusterPassFull_p60.sh</code>
    - and so on...

    - <code>sbatch mvt5_comp/51_ClusterPassFull_p15.sh</code>
    - <code>sbatch mvt5_comp/52_ClusterPassFull_p30.sh</code>
    - <code>sbatch mvt5_comp/53_ClusterPassFull_p60.sh</code>
    - and so on...  



6. **Generating Full Simulation Run Plots:** Once all full simulations are complete, you can generate visualizations of the full simulation runs. For this, follow all the indications of STEP 4 above, ensuring that you change the value of the variable <code>sim_type</code> in the files:

    - <code>Rscript 3_twosample_plots_sizepower_mvn_supps_Bcomp.r</code> 
    - <code>Rscript 4_twosample_plots_sizepower_mvn_supps_Norm.r</code> 
    - <code>Rscript 5_twosample_plots_sizepower_mvt5_comparison.r</code> 

    from <code>sim_type <- "exps"</code> to <code>sim_type <- "full"</code>. The resulting plots are saved in the <code>./results/</code> directory.

    **NOTE 1:** In all the figure-generating scripts described above, there is a variable <code>sim_type</code> which needs to be adjusted depending on whether the simulation is a reduced experiment or a full simulation run. Ensure you modify this value accordingly to obtain the appropriate results. 

    **NOTE 2:** The object <code>sim_type</code> is defined more than once within each of the plotting files. Before you run, ensure that you have set the value of <code>sim_type</code> to <code>sim_type <- "exps"</code> in all instances before running. 