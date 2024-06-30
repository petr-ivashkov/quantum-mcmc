#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --job-name="bo_cluster_TEST_100"
#SBATCH --mem-per-cpu=1024

export OMP_NUM_THREADS=48;

source $HOME/miniconda/bin/activate

conda activate qmcmc

# papermill ./bo_cluster/reverse_annealing_bo_n4.ipynb ./bo_cluster/reverse_annealing_bo_n4.ipynb
papermill ./bo_cluster/reverse_annealing_bo_n5.ipynb ./bo_cluster/reverse_annealing_bo_n5.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n6.ipynb ./bo_cluster/reverse_annealing_bo_n6.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n7.ipynb ./bo_cluster/reverse_annealing_bo_n7.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n8.ipynb ./bo_cluster/reverse_annealing_bo_n8.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n9.ipynb ./bo_cluster/reverse_annealing_bo_n9.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n10.ipynb ./bo_cluster/reverse_annealing_bo_n10.ipynb

# papermill ./standard_optim_cluster/reverse_annealing_cobyla_n4.ipynb ./standard_optim_cluster/reverse_annealing_cobyla_n4.ipynb


conda deactivate