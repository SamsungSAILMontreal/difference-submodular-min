#!/bin/sh

#SBATCH -o outputLogs/dsm_job-%A-%a.out
#SBATCH -a 0-125 
#SBATCH  --cpus-per-task 3
#SBATCH --time=1-00:00 # use 1 day for speech dataset, and 7 for mushroom dataset
#SBATCH --mem=32000M  

echo STARTING AT 'date'

echo "running on: "
hostname

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load matlab/2022a 

cd ~/diff-sub-min/Code

matlab -nodisplay -nodesktop -r "dsm_test(3, 'speech', 1, $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, 'results'); exit"
matlab -nodisplay -nodesktop -r "dsm_test(3, 'mushroom', 1e-4, $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, 'results'); exit"

echo "FINISHING AT `date`"
