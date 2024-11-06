#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --job-name="bo_n3"
#SBATCH --mem-per-cpu=1024

export OMP_NUM_THREADS=48;

source $HOME/miniconda/bin/activate
conda activate qmcmc

# papermill ./phase_transition/phase_diagram_quantum_SK_n4.ipynb ./phase_transition/phase_diagram_quantum_SK_n4.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n5.ipynb ./phase_transition/phase_diagram_quantum_SK_n5.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n6.ipynb ./phase_transition/phase_diagram_quantum_SK_n6.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n7.ipynb ./phase_transition/phase_diagram_quantum_SK_n7.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n8.ipynb ./phase_transition/phase_diagram_quantum_SK_n8.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n9.ipynb ./phase_transition/phase_diagram_quantum_SK_n9.ipynb

# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n3.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n3.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n4.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n4.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n5.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n5.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n6.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n6.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n7.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n7.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n8.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n8.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n9.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n9.ipynb

# papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n4.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n4.ipynb
# papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n5.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n5.ipynb
# papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n6.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n6.ipynb
# papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n7.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n7.ipynb
# papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n8.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n8.ipynb
# papermill ./grid_search_gamma_time/grid_search_time_gamma_parallel_n9.ipynb ./grid_search_gamma_time/grid_search_time_gamma_parallel_n9.ipynb

# papermill ./scale_factors/scale_factors_random.ipynb ./scale_factors/scale_factors_random.ipynb
# papermill ./scale_factors/scale_factors_local.ipynb ./scale_factors/scale_factors_local.ipynb
# papermill ./scale_factors/scale_factors_quantum.ipynb ./scale_factors/scale_factors_quantum.ipynb
# papermill ./scale_factors/scale_factors_fixed_gamma.ipynb ./scale_factors/scale_factors_fixed_gamma.ipynb
# papermill ./scale_factors/scale_factors_fixed_gamma_fixed_time.ipynb ./scale_factors/scale_factors_fixed_gamma_fixed_time.ipynb
papermill ./scale_factors/scale_factors_ra.ipynb ./scale_factors/scale_factors_ra.ipynb

# papermill ./schedule_optimization/reverse_annealing_bo_n4.ipynb ./schedule_optimization/reverse_annealing_bo_n4.ipynb
# papermill ./schedule_optimization/reverse_annealing_bo_n5.ipynb ./schedule_optimization/reverse_annealing_bo_n5.ipynb
# papermill ./schedule_optimization/reverse_annealing_bo_n6.ipynb ./schedule_optimization/reverse_annealing_bo_n6.ipynb
# papermill ./schedule_optimization/reverse_annealing_bo_n7.ipynb ./schedule_optimization/reverse_annealing_bo_n7.ipynb
# papermill ./schedule_optimization/reverse_annealing_bo_n8.ipynb ./schedule_optimization/reverse_annealing_bo_n8.ipynb
# papermill ./schedule_optimization/reverse_annealing_bo_n9.ipynb ./schedule_optimization/reverse_annealing_bo_n9.ipynb

# papermill ./time_dependency/optimal_t_vs_n.ipynb ./time_dependency/optimal_t_vs_n.ipynb

conda deactivate