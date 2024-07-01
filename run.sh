#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=10:00:00
#SBATCH --job-name="scale_factors"
#SBATCH --mem-per-cpu=4096

export OMP_NUM_THREADS=48;

source $HOME/miniconda/bin/activate

conda activate qmcmc

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


# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n4.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n4.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n5.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n5.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n6.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n6.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n7.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n7.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n8.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n8.ipynb
# papermill ./binder_cumulant/binder_cumulant_quantum_SK_n9.ipynb ./binder_cumulant/binder_cumulant_quantum_SK_n9.ipynb

papermill ./scale_factors/scale_factors_random.ipynb ./scale_factors/scale_factors_random.ipynb
papermill ./scale_factors/scale_factors_local.ipynb ./scale_factors/scale_factors_local.ipynb
papermill ./scale_factors/scale_factors_quantum.ipynb ./scale_factors/scale_factors_quantum.ipynb
papermill ./scale_factors/scale_factors_fixed_gamma.ipynb ./scale_factors/scale_factors_fixed_gamma.ipynb

conda deactivate