# Path
from src.path import projectdir

# Linear algebra
import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse.linalg as sparse_la
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.interpolate import interp1d

# Plotting and output
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
import seaborn as sns

# Useful stuff
import time
import json
import pickle
from tqdm import tqdm
import joblib

# Define colors with names
colors = {
    "red": "#e41a1c",
    "blue": "#377eb8",
    "green": "#4daf4a",
    "purple": "#984ea3",
    "orange": "#ff7f00",
    "yellow": "#ffff33",
    "brown": "#a65628",
    "pink": "#f781bf",
    "grey": "#999999",
    "black": "#000000",
    "white": "#ffffff"
}

plt.rcParams.update(
    {
        "xtick.direction": "in",
        "ytick.direction": "out",
        "ytick.right": False,
        "xtick.top": False,
        "ytick.left": True,
        "xtick.bottom": False,
        "figure.facecolor": "1",
        "savefig.facecolor": "1",
        "savefig.dpi": 600,
        "figure.dpi": 600,
        "savefig.bbox": "tight",
        "font.size": 7,
        "font.family": "serif",
        "lines.markersize": 6,
        "lines.linewidth": 1,
        'axes.axisbelow' : True
    }
)

# Golden ratio
figure_size_x = 6.0462
figure_size_y = figure_size_x/1.618

# Conversion funtions
def int_to_bin(i, n):
    '''Convert a given integer to a bitstring of fixed length using 0 -> 0..00 convention.'''
    return bin(i)[2:].zfill(n)

def int_to_bin_arr(i, n):
    '''Convert a given integer to a bitstring of fixed length using 0 -> 0..00 convention.'''
    return np.array([int(b) for b in int_to_bin(i, n)])

def bin_to_int(s):
    '''Convert a given bitstring to an integer.'''
    return int(s,2)

def bin_to_spin(b):
    '''Convert a binary string to a spin configuration using 0 -> +1 convention.'''
    if len(b) == 1: return 1-2*int(b)
    else:
        s = [1-2*int(c) for c in b]
        return s

def int_to_spin(i, n):
    '''Convert an integer to a spin configuration of specified length'''
    b = int_to_bin(i, n)
    return bin_to_spin(b)

def get_computational_basis_state(i, n):
    '''Return the i-th computational basis state vector of dimension 2^n.'''
    state = np.zeros(2**n, dtype=complex)
    state[i] = 1
    return state

# Problem instance
class IsingModel:
    '''
    Represents an Ising model with given parameters J and h.

    Attributes:
    - n (int): Number of spins.
    - h (numpy.ndarray): External field.
    - J (numpy.ndarray): Interaction matrix.
    - T (float): Temperature (default 1).
    - J_rescaled (numpy.ndarray): Rescaled interaction matrix.
    - h_rescaled (numpy.ndarray): Rescaled external field.
    - E (numpy.ndarray): Energies of all possible spin configurations.
    '''
    def __init__(self, J, h):
        self.h = h
        self.n = h.size
        self.J = J

        assert np.allclose(J.T, J), 'J must be symmetric.'
        assert np.all(np.diag(J) == 0), 'Diagonal elements of J must be zero.'

        self.alpha = np.sqrt(self.n / (0.5*la.norm(J, ord='fro')**2 + la.norm(h, ord=2)**2))
        self.J_rescaled = self.J * self.alpha 
        self.h_rescaled = self.h * self.alpha

        self.E = np.zeros((2**self.n))
        for i in range(2**self.n):
            s = int_to_spin(i, self.n)
            self.E[i] = -0.5*(s @ self.J @ s) - s @ self.h
        self.gs = int_to_bin(np.argmin(self.E), self.n)
        self.gs_deg = np.count_nonzero(np.isclose(self.E, min(self.E)))
        self.E_rescaled = self.E.copy() * self.alpha

    @classmethod
    def from_coefficients(cls, n, coefficients):
        '''
        Instantiates an Ising model from an array of coefficients.
        
        Parameters:
        - coefficients (numpy.ndarray): 1D array of n + n(n-1)/2 coefficients, where the first n elements
                                        correspond to h (external fields) and the remaining correspond to
                                        the upper triangular elements of the symmetric interaction matrix J.
        '''
        #  Check that the length of the coefficients array is consistent
        assert len(coefficients) == n+(n*(n-1))//2, \
            f"Expected {n+(n*(n-1))//2} coefficients, but got {len(coefficients)}."

        h = np.asarray(coefficients[:n])
        J_upper_tri = coefficients[n:]
        J = np.zeros((n, n))

        # Fill the upper triangular part of J with the interaction coefficients
        upper_tri_indices = np.triu_indices(n, k=1)
        J[upper_tri_indices] = J_upper_tri
        J += J.T
        return cls(J, h)

    def __str__(self):
        info = f"Ising model information:\n"
        info += f"Number of spins: {self.n}\n"
        info += f"External fields (h): {self.h}\n"
        info += f"Interaction matrix (J):\n{self.J}\n"
        info += f"Alpha: {self.alpha}\n"
        info += f"Ground state energy: {min(self.E)}\n"
        info += f"Ground state: |{self.gs}> = |{bin_to_int(self.gs)}>\n"
        info += f"Ground state degeneracy: {self.gs_deg}\n"
        info += f"Energy spectrum: {self.E}\n"
        return info

class RandomIsingModel(IsingModel):
    '''
    Represents a randomly generated Ising model with a given number of qubits.
    '''
    def __init__(self, n, local_fields=True, seed=None):
        if local_fields:
            h = generate_random_h(n, seed=seed)
        else: 
            h = np.zeros(n)
        J = generate_random_J(n, seed=seed)
        super().__init__(J, h)

# Useful quantum definitions - taken from D. Layden et al [Nature 619, 282–287 (2023)] 
def my_kron(arr_list):
    '''Takes a list of numpy arrays [A1, A2, ...] and computes their tensor product A1 (x) A2 (x) ...'''
    if len(arr_list) == 1:
        return arr_list[0]
    else:
        return np.kron(arr_list[0], my_kron(arr_list[1:]))

def X(i, n):
    '''n-qubit Pauli X. Acts as sigma_x on qubit i and as the identity on the rest.'''
    if i<0 or i>=n or n<1:
        raise ValueError('Bad value of i and/or n.')
    X_list = [np.array([[0,1],[1,0]]) if j==i else np.eye(2) for j in range(n)]
    return my_kron(X_list)

def Y(i, n):
    '''n-qubit Pauli Y. Acts as sigma_y on qubit i and as the identity on the rest.'''
    if i<0 or i>=n or n<1:
        raise ValueError('Bad value of i and/or n.')
    Y_list = [np.array([[0,-1j],[1j,0]]) if j==i else np.eye(2) for j in range(n)]
    return my_kron(Y_list)

def Z(i, n):
    '''n-qubit Pauli Z. Acts as sigma_z on qubit i and as the identity on the rest.'''
    if i<0 or i>=n or n<1:
        raise ValueError('Bad value of i and/or n.')
    Z_list = [np.array([[1,0],[0,-1]]) if j==i else np.eye(2) for j in range(n)]
    return my_kron(Z_list)

# Functions to generate a random instance 
def generate_random_J(n, seed=None):
    '''Generate a random symmetric coupling matrix for n sites, sampled from N(0,1).'''
    np.random.seed(seed)
    J = np.triu(np.random.standard_normal((n,n)), k=1) # No self-interactions
    J = J + J.T # Ensure symmetry
    return J

def generate_random_h(n, seed=None):
    '''Generate random external fields for n sites, sampled from N(0,1).'''
    if seed is not None:
        seed += int((n+1)*n/2)
    np.random.seed(seed) 
    return np.random.standard_normal(n)

# Functions for generating a classical proposal
def hamming(s1, s2):
    '''Calculate the Hamming distance between two bit strings.'''
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def get_proposal_mat_random(m):
    '''Get a random proposal matrix for a given Ising model. Every configuration is proposed with equal probability.'''
    return np.ones((2**m.n, 2**m.n)) / 2**m.n

def get_proposal_mat_local(m, hamming_radius=1):
    '''Get a local proposal matrix for a given Ising model. Only configurations
       within a specified Hamming radius can be proposed.'''
    assert hamming_radius <= m.n
    # Note: diagonal elements of proposal_mat must be zero
    d= 2**m.n
    p_propose = 1 / sum([scipy.special.binom(m.n,r) for r in range(1,hamming_radius+1)]) # proposal probability
    proposal_mat = np.zeros((d,d))
    for i in range(d):
        s_i = int_to_bin(i,m.n)
        for j in range(d):
            s_j = int_to_bin(j,m.n)
            if hamming(s_i, s_j) > hamming_radius: continue # only local spinflips can be proposed
            if i == j: continue
            proposal_mat[i,j] = p_propose
    assert is_stochastic(proposal_mat), 'Something went wrong: Matrix is not stochastic.' 
    return proposal_mat

# Functions for generating a quantum proposal
def get_basis_state(i,n):
    '''Get the computational basis state |i> of n qubits.'''
    assert (i < 2**n) & (i >= 0), 'Bad state.'
    psi = np.zeros(1<<n)
    psi[i] = 1
    return psi

def get_mixing_hamiltonian(n):
    '''Simple sparse mixing Hamiltonian.'''
    return sum([X(i,n) for i in range(n)])

# Precompute mixing Hamiltonians 
H_mixer_list = [get_mixing_hamiltonian(n) for n in range(1,11)]

def get_proposal_mat_quantum(m, gamma=0.7, t=1):
    '''Get a quantum proposal matrix for a given Ising model.'''
    H_ising = np.diag(m.E_rescaled)
    H_mixer = H_mixer_list[m.n-1]
    H = (1-gamma)*H_ising + gamma*H_mixer
    U = sparse_la.expm(-1j * H * t) # time evoluiton operator
    proposal_mat = np.abs(U)**2
    return proposal_mat

def get_proposal_mat_quantum_avg(m, gamma_lims=(0.25, 0.6), gamma_steps=20, time_lims=(2,12)):
    '''Get a quantum proposal matrix for a given Ising model averaged over time and gamma.'''
    # Constructing the Ising and mixing Hamiltonians
    H_ising = np.diag(m.E_rescaled)
    H_mixer = H_mixer_list[m.n-1]

    t_min = time_lims[0]
    t_max = time_lims[1]

    # approximate integral over gamma by Riemann sum with gamma_steps (D. Layden et al [Nature 619, 282–287 (2023)])
    gamma_starts, step_size = np.linspace(gamma_lims[0], gamma_lims[1], num=gamma_steps, endpoint=False, retstep=True)
    gamma_range = gamma_starts + step_size/2

    d = 2**m.n
    proposal_mat_arr = []
    for gamma in gamma_range:
        H = (1-gamma)*H_ising + gamma*H_mixer
        vals, vecs = la.eigh(H)

        # Construct transition weight matrix for to time averaging (analytical expression)
        transition_weight_mat = np.zeros((d,d))
        for i in range(d):
            for j in range(i+1,d):
                dlambda = vals[i] - vals[j]
                if dlambda != 0:
                    transition_weight_mat[i,j] = (np.sin(dlambda*t_max) - np.sin(dlambda*t_min) ) / (dlambda * (t_max-t_min))
                else: transition_weight_mat[i,j] = 1
        transition_weight_mat = transition_weight_mat + transition_weight_mat.T + np.eye(d)
        # Construct the time-averaged proposal matrix
        Q = np.zeros((d,d))
        for i in range(d):
            for j in range(i,d):
                A = np.outer(vecs[i], vecs[i])
                B = np.outer(vecs[j], vecs[j])
                Q[i,j] = np.sum(A*B*transition_weight_mat)
        Q = Q + Q.T - np.diag(np.diag(Q))
        proposal_mat_arr.append(Q)

    # Take an average over gammas
    proposal_mat = np.mean(proposal_mat_arr, axis=0) 
    assert is_stochastic(proposal_mat), 'Proposal matrix is not stochastic.'
    return proposal_mat

def get_proposal_mat_quantum_layden(m, gamma_lims=(0.25, 0.6), gamma_steps=20, time_lims=(2,20)):
    '''
    Returns a 2**n * 2**n stochastic proposal matrix for our quantum method of suggesting moves,
    with no Trotter, gate or SPAM errors.

    Taken from D. Layden et al [Nature 619, 282–287 (2023)].

    This function is optimized and is significantly faster than my <get_proposal_mat_quantum_avg>.
    ''' 
    def cont_eig(Dlambda):
        t_0, t_f = time_lims[0], time_lims[1] # evolution time t ~ unif(t_0, t_f)
        x = np.sin(Dlambda*t_f) - np.sin(Dlambda*t_0) # from analytical expression for transition probabilities
        return np.divide(2*x/(t_f-t_0), Dlambda, out=np.ones_like(Dlambda), where=(Dlambda!=0) )
    
    n = m.n 

    H_ising = np.diag(m.E_rescaled)
    H_mixer = H_mixer_list[m.n-1]

    d = 2**n
    a = np.arange(d**2)
    mask = (a//d >= a%d)
    ones = np.ones(d)

    # approximate integral over gamma by Riemann sum with gamma_steps (D. Layden et al [Nature 619, 282–287 (2023)])
    gamma_starts, step_size = np.linspace(gamma_lims[0], gamma_lims[1], num=gamma_steps, endpoint=False, retstep=True)
    gamma_range = gamma_starts + step_size/2

    prop_list = []
    for gamma in gamma_range:
        H = (1-gamma)*H_ising + gamma*H_mixer
        vals, vecs = la.eigh(H)
        vals_diff = (np.kron(vals, ones) - np.kron(ones, vals))[mask]
        M = la.khatri_rao(vecs.T, vecs.T)[mask]
        Q = M.T * cont_eig(vals_diff) @ M 
        prop_list.append(Q)
    
    proposal_mat = sum(prop_list)/gamma_steps
    return proposal_mat

# MCMC functions
def is_stochastic(P):
    '''Checks if matrix P is left stochastic, i.e., columns add up to one.'''
    col_sums = np.sum(P, axis=0)
    return np.all(P >= -1e-8) and np.allclose(col_sums, 1)

def get_transition_matrix_old_version(m, T, proposal_mat):
    '''
    Generate a MC transition probability matrix using Metropolis-Hastings algorithm. 
    This is the first version of this function. Use <get_transition_matrix> for better performance.
    '''
    d = 2**m.n
    P = np.zeros((d,d))
    for i in range(d):
        for j in range(d):
            if i == j: continue # sum for j != i
            dE = m.E[i] - m.E[j] # E_new - E_old
            if dE > 0: A = np.exp(-dE / T) # MH acceptance
            else: A = 1 # only compute the exponential if dE > 0 to avoid overflows
            P[i,j] = A * proposal_mat[i,j]
    for i in range(d):
        P[i,i] = 1 - sum(P[:,i]) # sum for j == i for normalization
    assert is_stochastic(P), 'Something went wrong: P is not stochastic.'
    return P

def get_transition_matrix(m, T, proposal_mat):
    '''Generate a MC transition probability matrix using Metropolis-Hastings algorithm.'''
    d = 2**m.n
    E_tiled = np.tile(m.E, (d,1))
    dE = E_tiled.T - E_tiled # dE[new,old] = E_new - E_old
    if T != 0:
        A = np.exp(-dE / T, where=(dE > 0), out=np.ones_like(dE)) # only compute the exponential if dE > 0 to avoid overflows
    else: 
        A = np.where(dE > 0, 0., 1.) # reject all energy-increasing steps
    P = A * proposal_mat # MH step

    np.fill_diagonal(P, 0)
    diag = np.ones(d) - np.sum(P, axis=0)
    P = P + np.diag(diag)
    return P

def get_delta(P):
    '''Calculate the spectral gap of the transition matrix P.'''
    vals = np.abs(la.eigvals(P)) # sparse_la.eigs has problems converging for n > 5
    sorted_vals = sorted(vals, reverse=True)
    assert np.isclose(sorted_vals[0], 1), f'The largest eigenvalue {sorted_vals[0]} differs from 1.'
    return sorted_vals[0] - sorted_vals[1]

def get_transition(P, state):
    '''Generate a random transition from a given state using a given transition matrix.'''
    if type(state) == type('s'): state = bin_to_int(state)
    return np.random.choice(P.shape[0], p=P[:,state])

def get_trajectory(P, num_moves, initial_state=None): 
    '''Generate a trajectory of <num_moves> moves based on the transition matrix <P>.'''
    if initial_state is None: initial_state = np.random.randint(P.shape[0])
    if type(initial_state) == type('s'): initial_state = bin_to_int(initial_state)
    trajectory = np.empty(num_moves+1, dtype='int')
    trajectory[0] = initial_state
    for i in range(1, num_moves+1):
        trajectory[i] = get_transition(P, trajectory[i-1])
    return trajectory

def samples_to_counts(samples, dtype='str'):
    'Convert sample list to a dictionary of counts'
    counts = {}
    if dtype == 'str':
        for sample in samples:
            if sample in counts: counts[sample] += 1
            else: counts[sample] = 1
    if dtype == 'int':
        for sample in samples:
            if bin_to_int(sample) in counts: counts[bin_to_int(sample)] += 1
            else: counts[bin_to_int(sample)] = 1  
    return counts

# TDSE solvers
def state_tdse(t, psi, hamiltonian):
    '''
    Computes the time derivative of the state vector psi 
    according to the time-dependent Schrödinger equation.
    '''
    return -1j * hamiltonian(t) @ psi

def operator_tdse(t, U, hamiltonian):
    '''
    Computes the time derivative of the unitary operator U 
    according to the time-dependent Schrödinger equation.
    '''
    return -1j * hamiltonian(t) @ U

def solve_state_tdse(hamiltonian, t_max, psi_0, return_intermediate=False, t_eval=None, method='RK45', rtol=1e-10):
    '''
    Solves the time-dependent Schrödinger equation for the 
    state vector psi over the interval [0, t_max].
    '''
    sol = solve_ivp(state_tdse,
                    [0, t_max],
                    y0=psi_0,
                    method=method,
                    t_eval=t_eval,
                    rtol=rtol,
                    args=(hamiltonian,) # args for state_tdse()
                    )
    psi_t = sol.y.T # shape of sol.y is d x n_val_points
    if return_intermediate:
        return psi_t
    else: 
        return psi_t[-1]
    
def solve_operator_tdse(hamiltonian, t_max, return_intermediate=False, t_eval=None, method='RK45', rtol=1e-10):
    '''
    Solves the time-dependent Schrödinger equation for a 
    unitary operator over the interval [0, t_max].
    '''
    def flattened_operator_tdse(t, U_flat, hamiltonian):
        # Returns a flattened U for solving TDSE using solve_ivp()
        return operator_tdse(t, U_flat.reshape(d, d), hamiltonian).flatten()
    
    d = len(hamiltonian(0))
    U0 = np.eye(d, dtype=complex).flatten()
    sol = solve_ivp(flattened_operator_tdse,
                    [0, t_max],
                    y0=U0,
                    method=method,
                    t_eval=t_eval,
                    rtol=rtol,
                    args=(hamiltonian,) # args for operator_tdse()
                    )
    U_t = sol.y.T.reshape((-1,d,d)) # shape of sol.y is (d x d) x n_val_points
    if return_intermediate:
        return U_t
    else: 
        return U_t[-1]


# Annealing related functions

RAMP_UP_TIME_FRACTION = 0.05

def get_symmetric_schedule(ramp_up_schedule, start_at_zero=True):
    '''Mirror and concatenate the given ramp up schedule to create a symmetric schedule.'''
    if start_at_zero:
        ramp_up_schedule = np.concatenate(([0], ramp_up_schedule))
    return np.concatenate([ramp_up_schedule, ramp_up_schedule[::-1][1:]])

def get_proposal_mat_ra(m, schedule_interpolator, t_max, assert_symmetry=True):
    '''
    Computes the proposal matrix using the given schedule function 
    for the reverse annealing process with the total annealing time t_max.

    - schedule() should be a symmetric function on the interval [0,1].
    '''
    H_ising = np.diag(m.E_rescaled)
    H_mixer = H_mixer_list[m.n-1]
    def s(t): return t/t_max # annealing fraction
    def hamiltonian(t): return (1-schedule_interpolator(s(t)))*H_ising + schedule_interpolator(s(t))*H_mixer
    U_t = solve_operator_tdse(hamiltonian, t_max)
    if assert_symmetry: 
         assert np.allclose(U_t, U_t.T), 'Proposal matrix is not symmetric.'
    return np.abs(U_t)**2

def get_schedule_interpolator(schedule, kind='linear'):
    '''Creates an interpolating function based on the schedule.'''
    n_points = len(schedule)
    # s = np.linspace(RAMP_UP_TIME_FRACTION,1-RAMP_UP_TIME_FRACTION, n_points-2)
    # s = np.concatenate(([0],s,[1]))
    s = np.linspace(0,1,n_points)
    return interp1d(s, schedule, kind=kind)

# Order parameters for quantum SK model
def magnetization(H, T):
    d = len(H)
    n = int(np.log2(d))
    mag_basis = np.zeros(d)
    for i in range(d):
        s = int_to_bin(i, n)
        mag_basis[i] = (s.count('0') - s.count('1')) / n
    E, vecs = la.eigh(H)
    if T == 0:
        boltzmann_factors = np.array(E == min(E), dtype=int) / sum(np.array(E == min(E), dtype=int))
    else:
        boltzmann_factors = scipy.special.softmax(-E / T)
    mag = 0
    for i, vec in enumerate(vecs.T):
        mag += mag_basis @ vec**2 * boltzmann_factors[i]
    return mag

def qea(H, T):
    '''Edward-Anderson parameter'''
    d = len(H)
    n = int(np.log2(d))
    E, vecs = la.eigh(H)
    if T == 0:
        boltzmann_factors = np.array(E == min(E), dtype=int) / sum(np.array(E == min(E), dtype=int))
    else:
        boltzmann_factors = scipy.special.softmax(-E / T)

    qea_val = 0
    for i in range(n):
        sigma_i = Z(i, n)
        s = 0
        for vec_index, vec in enumerate(vecs.T):
            s += vec @ sigma_i @ vec.T * boltzmann_factors[vec_index]
        qea_val += s**2

    return qea_val / n

# Plotting functions
def plot_schedule(schedule, schedule_interpolator = None):
    '''Plots annealing schedule with optional interpolation.'''
    n_points = len(schedule)
    # s = np.concatenate(([0],np.linspace(RAMP_UP_TIME_FRACTION,1-RAMP_UP_TIME_FRACTION, n_points-2),[1]))
    s = np.linspace(0,1,n_points)
    if schedule_interpolator is None:
        plt.plot(s, schedule, '.--', color = colors['red'])
    else:
        x_cont = np.linspace(0,1,100)
        plt.plot(x_cont, schedule_interpolator(x_cont), '--', color = colors['grey'])
        plt.plot(s, schedule, '.', color = colors['red'])
    plt.xlabel('s')
    plt.ylabel('$\gamma(s)$')
    plt.title('$H(\gamma) = (1-\gamma) H_{prob} + \gamma H_{mix}$')
    plt.show()

def display_video(video, fps=30, cmap='coolwarm'):
    """
    Display a 3D numpy array as a video in a Jupyter notebook.

    Parameters:
    array (numpy.ndarray): A 3D numpy array where each slice along the 3rd axis is a frame.
    fps (int): Frames per second, controls the speed of the video.
    """
    n_frames = video.shape[0]
    fig, ax = plt.subplots()
    img = ax.imshow(video[0,:,:], 
                    interpolation='nearest', 
                    origin='lower',
                    cmap=cmap)
    for i in range(n_frames):
        frame = video[i,:,:]
        img.set_data(frame)
        display(fig)
        clear_output(wait=True)
        time.sleep(1 / fps)
    plt.close(fig)

# Other useful functions
def save_in_json(data, file_path):
    '''Save the dictionary as a JSON file.'''
    with open(projectdir+file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)
    # Check if writing was successful
    with open(projectdir+file_path, 'r') as json_file:
        data_loaded = json.load(json_file)
    assert data_loaded == data, 'An error occured when saving JSON.'

def load_from_json(file_path):
    '''Load data from a JSON file and return it as a dictionary.'''
    with open(projectdir+file_path, 'r') as json_file:
        data = json.load(json_file)
    return data