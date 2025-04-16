#!/bin/sh
#SBATCH --job-name=MVN960
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=10-00:00:00
#SBATCH --output=mvn_comp/output_full/output_p960_%a.out
#SBATCH --array=100-6099

pwd; hostname; date

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo This is task $((SLURM_ARRAY_TASK_ID)),
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo 1 Importing R and Rscript:
module load r/4.4.0
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo 3 Running mvn_comp/31_SimulationScript.R with input $((SLURM_ARRAY_TASK_ID + 36000)) 2
Rscript mvn_comp/31_SimulationScript.R $((SLURM_ARRAY_TASK_ID + 36000)) 2
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

date

