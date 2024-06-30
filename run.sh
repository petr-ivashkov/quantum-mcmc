#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --job-name="bo_cluster_TEST"
#SBATCH --mem-per-cpu=1024

export OMP_NUM_THREADS=48;

source $HOME/miniconda/bin/activate

conda activate qmcmc

# papermill ./bo_cluster/reverse_annealing_bo_n4.ipynb ./bo_cluster/reverse_annealing_bo_n4.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n5.ipynb ./bo_cluster/reverse_annealing_bo_n5.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n6.ipynb ./bo_cluster/reverse_annealing_bo_n6.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n7.ipynb ./bo_cluster/reverse_annealing_bo_n7.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n8.ipynb ./bo_cluster/reverse_annealing_bo_n8.ipynb
papermill ./bo_cluster/reverse_annealing_bo_n9.ipynb ./bo_cluster/reverse_annealing_bo_n9.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n10.ipynb ./bo_cluster/reverse_annealing_bo_n10.ipynb

conda deactivate