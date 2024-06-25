#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=20:00:00
#SBATCH --job-name="grid_search_time_gamma_parallel_n5-10"
#SBATCH --mem-per-cpu=1024

export OMP_NUM_THREADS=48;

source /cluster/home/pivashkov/virtualenvs/qmcmc/bin/activate

papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n5.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n5.ipynb
papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n6.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n6.ipynb
papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n7.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n7.ipynb
papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n8.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n8.ipynb
papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n9.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n9.ipynb
papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n10.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n10.ipynb

deactivate