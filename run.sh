#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
#SBATCH --job-name="time"
#SBATCH --mem-per-cpu=2056

export OMP_NUM_THREADS=48;

# papermill ./bo_cluster/reverse_annealing_bo_n4.ipynb ./bo_cluster/reverse_annealing_bo_n4.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n5.ipynb ./bo_cluster/reverse_annealing_bo_n5.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n6.ipynb ./bo_cluster/reverse_annealing_bo_n6.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n7.ipynb ./bo_cluster/reverse_annealing_bo_n7.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n8.ipynb ./bo_cluster/reverse_annealing_bo_n8.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n9.ipynb ./bo_cluster/reverse_annealing_bo_n9.ipynb
# papermill ./bo_cluster/reverse_annealing_bo_n10.ipynb ./bo_cluster/reverse_annealing_bo_n10.ipynb

# papermill ./phase_transition/phase_diagram_quantum_SK_n4.ipynb ./phase_transition/phase_diagram_quantum_SK_n4.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n5.ipynb ./phase_transition/phase_diagram_quantum_SK_n5.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n6.ipynb ./phase_transition/phase_diagram_quantum_SK_n6.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n7.ipynb ./phase_transition/phase_diagram_quantum_SK_n7.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n8.ipynb ./phase_transition/phase_diagram_quantum_SK_n8.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n9.ipynb ./phase_transition/phase_diagram_quantum_SK_n9.ipynb
# papermill ./phase_transition/phase_diagram_quantum_SK_n10.ipynb ./phase_transition/phase_diagram_quantum_SK_n10.ipynb

# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n3.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n3.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n4.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n4.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n5.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n5.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n6.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n6.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n7.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n7.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n8.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n8.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n9.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n9.ipynb

# papermill ./scale_factors/scale_factors_random.ipynb ./scale_factors/scale_factors_random.ipynb
# papermill ./scale_factors/scale_factors_local.ipynb ./scale_factors/scale_factors_local.ipynb
# papermill ./scale_factors/scale_factors_quantum.ipynb ./scale_factors/scale_factors_quantum.ipynb
# papermill ./scale_factors/scale_factors_fixed_gamma.ipynb ./scale_factors/scale_factors_fixed_gamma.ipynb
# papermill ./scale_factors/scale_factors_ra.ipynb ./scale_factors/scale_factors_ra.ipynb
papermill ./time_dependency/optimal_t_vs_n.ipynb ./time_dependency/optimal_t_vs_n.ipynb
