#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=20:00:00
#SBATCH --job-name="phase_diagram_quantum_SK_n5_10"
#SBATCH --mem-per-cpu=1024

export OMP_NUM_THREADS=48;

source /cluster/home/pivashkov/virtualenvs/qmcmc/bin/activate

papermill ./phase_transition/phase_diagram_quantum_SK_n5.ipynb ./phase_transition/phase_diagram_quantum_SK_n5.ipynb
papermill ./phase_transition/phase_diagram_quantum_SK_n6.ipynb ./phase_transition/phase_diagram_quantum_SK_n6.ipynb
papermill ./phase_transition/phase_diagram_quantum_SK_n7.ipynb ./phase_transition/phase_diagram_quantum_SK_n7.ipynb
papermill ./phase_transition/phase_diagram_quantum_SK_n8.ipynb ./phase_transition/phase_diagram_quantum_SK_n8.ipynb
papermill ./phase_transition/phase_diagram_quantum_SK_n9.ipynb ./phase_transition/phase_diagram_quantum_SK_n9.ipynb
papermill ./phase_transition/phase_diagram_quantum_SK_n10.ipynb ./phase_transition/phase_diagram_quantum_SK_n10.ipynb

deactivate