#!/bin/bash

#SBATCH -A es_pelli
#SBATCH -J casCHades
#SBATCH --time=00-4:00:00
#SBATCH --output=log/run_%A/simul_%A.out
#SBATCH --error=log/run_%A/simul_%A.error
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=0-999  # Adjust the array size accordingly

module load gcc/9.3.0 python/3.11.2

# Calculate the total number of simulations
total_simulations=4000 # n simulations per task, starting from ID 0

# Calculate the number of simulations per task
simulations_per_task=4

# Calculate the start ID for the current job array task
start_id=$((SLURM_ARRAY_TASK_ID * simulations_per_task))
end_id=$((start_id + simulations_per_task - 1))

# Loop through each task
for ((sim_id=start_id; sim_id<=end_id; sim_id++))
do
    python Combined.py -i $start_id
done

