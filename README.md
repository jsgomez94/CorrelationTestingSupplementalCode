# CorrTestingCodeSubmission
This repository contains all the proof-read and commented code to recreate simulations, figures and tables of the paper "Effective Permutation Tests for Differences Across Multiple High-Dimensional Correlation Matrices."

Here, we simply describe the overall structure of the repository. To fully reproduce the results, further instructions are necessary. For this end, we provide further README files when necessary. The structure is the following:

- <code>./real_data/</code>: directory containing the code and data for performing the real data analysis we describe in Section 6 of our main paper. We provide a README file providing further instructions on how to replicate the real data analysis. ***NOTE***: We do not provide the code and data for reproducing the analysis of brain activation measurements for Alzheimer's patients studied in Section 5.2 of our paper, due to privacy concerns regarding our dataset.

- <code>./simulations/</code>: directory containing the code for reproducing the simulation results of Section 4 in our main paper, and Sections S2, S3, S4, S5, S6, S7 and S8  of our Supplementary Materials. In order to replicate all simulation experiments and plots, you must run the simulation code on a linux cluster with SLURM scheduling system. We provide further description of how to reproduce the results in the README of this directory.S