# Quantum MCMC for Ising Spin Glasses

This repository is the official codebase accompanying [Phys. Rev. A **111**, 042615 (2025)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.111.042615). It explores quantum-inspired Markov chain Monte Carlo (MCMC) sampling for Ising spin-glass instances and performs a comprehensive analysis of the quantum-enhanced Monte Carlo method  [Nature **619**, 282–287 (2023)](https://www.nature.com/articles/s41586-023-06095-4).

## Repository Layout
- `src/functions.py` - Core library for building Ising instances, quantum operators, proposal matrices, Metropolis-Hastings transitions, time-dependent Schrodinger solvers, annealing schedules, observables, plotting helpers, and JSON I/O.
- `src/path.py` - Project root path helper (`projectdir`); update this to your local path if you persist outputs.
- Experiment notebooks live under `plots/`, `phase_transition/`, `binder_cumulant/`, `scale_factors/`, `grid_search_gamma_time/`, `schedule_optimization/`, `qaoa_data/`, `time_dependency/`, and `miscellaneous/`, covering phase diagrams, spectral gaps, Binder cumulants, schedule optimization, reverse annealing, QAOA baselines, and related analyses.
- `data/` - Precomputed JSON/Pickle datasets (grid searches, spectral gaps, Binder cumulants, temperature sweeps, scale factors, qubit sweeps).
- `figures/` - Exported plots (PNG/PDF) corresponding to the notebooks.
- `run.sh` - Example Slurm + `papermill` batch runner for reproducing notebooks on an HPC cluster.
- `environment.yml`, `requirements.txt` - Conda/Pip environments.

## Core Library (`src/functions.py`)
- **Bitstring/spin utilities**: `int_to_bin`, `int_to_bin_arr`, `bin_to_int`, `bin_to_spin`, `int_to_spin`, `get_computational_basis_state`.
- **Ising instances**: `IsingModel` (energy spectrum, rescaling, ground state), `RandomIsingModel`, and `from_coefficients` constructor; random generators `generate_random_J`, `generate_random_h`.
- **Pauli operators & mixing Hamiltonian**: Dense `X(i,n)`, `Y(i,n)`, `Z(i,n)` with recursive `my_kron`; `get_mixing_hamiltonian` and precomputed `H_mixer_list`.
- **Proposal matrices**:
  - Classical: `get_proposal_mat_random` (global), `get_proposal_mat_local` (Hamming ball), plus `hamming`.
  - Quantum: `get_proposal_mat_quantum` (time-evolved mixer/problem Hamiltonian blend), `get_proposal_mat_quantum_avg` (gamma/time averaging), `get_proposal_mat_quantum_layden` (vectorized/optimized implementation of the proposal scheme; follows the of Nature 619, 282–287 (2023)).
- **Metropolis-Hastings transitions**: `get_transition_matrix` (vectorized MH), `get_transition_matrix_old_version`, `is_stochastic`, `get_delta` (spectral gap), `get_transition`, `get_trajectory`, `samples_to_counts`.
- **TDSE solvers**: `solve_state_tdse`, `solve_operator_tdse` with helpers (`state_tdse`, `operator_tdse`); used for quantum proposal evolution and annealing schedules.
- **Annealing schedules**: `get_symmetric_schedule`, `get_schedule_interpolator`, `get_proposal_mat_ra` (reverse-annealing proposal from a schedule).
- **Order parameters**: `magnetization` and `qea` (Edwards-Anderson) for quantum Sherrington-Kirkpatrick analyses.
- **Plotting/visualization**: Matplotlib/Seaborn defaults, `plot_schedule`, and `display_video` for notebook-friendly animations.
- **Persistence**: `save_in_json`, `load_from_json` (uses `src/path.py`'s `projectdir`).

## Installation
```bash
conda env create -f environment.yml
conda activate qmcmc
```

## Typical Usage
```python
from src.functions import (
    RandomIsingModel, get_proposal_mat_quantum,
    get_transition_matrix, get_trajectory, get_delta
)

m = RandomIsingModel(n=5, seed=0)
proposal = get_proposal_mat_quantum(m, gamma_steps=20)
P = get_transition_matrix(m, T=0.5, proposal_mat=proposal)
gap = get_delta(P)
traj = get_trajectory(P, num_moves=1000)
```
For reverse annealing, build a schedule with `get_schedule_interpolator` and pass it to `get_proposal_mat_ra`.

## Citation
If you use this code or reproduce the proposal mechanisms, please cite the paper above (Phys. Rev. A **111**, 042615 (2025)).

## Notes
- Update `src/path.py` (`projectdir`) to your local clone if you save JSON outputs.