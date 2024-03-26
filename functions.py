import numpy as np
import time
import scipy
import scipy.linalg as la
import scipy.sparse.linalg as sparse_la

from IPython.display import display, clear_output
import matplotlib.pyplot as plt

# Conversion funtions
def int_to_bin(i, n):
    '''Convert a given integer to a bitstring of fixed length using 0 -> 0..00 convention.'''
    return bin(i)[2:].zfill(n)

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

        self.alpha = np.sqrt(self.n / (la.norm(J, ord='fro')**2 + la.norm(h, ord=2)**2))
        self.J_rescaled = self.J * self.alpha 
        self.h_rescaled = self.h * self.alpha

        self.E = np.zeros((2**self.n))
        for i in range(2**self.n):
            s = int_to_spin(i, self.n)
            self.E[i] = -0.5*(s @ self.J @ s) - s @ self.h
        self.E_rescaled = self.E.copy() * self.alpha

class RandomIsingModel(IsingModel):
    '''
    Represents a randomly generated Ising model with a given number of qubits.
    '''
    def __init__(self, n):
        h = generate_random_h(n)
        J = generate_random_J(n)
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
def generate_random_J(n):
    '''Generate a random symmetric coupling matrix for n sites, sampled from N(0,1).'''
    J = np.triu(np.random.standard_normal((n,n)), k=1) # No self-interactions
    J = J + J.T # Ensure symmetry
    return J

def generate_random_h(n):
    '''Generate random external fields for n sites, sampled from N(0,1).'''
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

def get_proposal_mat_quantum(m, gamma=0.7, t=1):
    '''Get a quantum proposal matrix for a given Ising model.'''
    H_zz = sum([-m.J_rescaled[i,j]*Z(i,m.n) @ Z(j,m.n) for i in range(m.n) for j in range(i+1,m.n)]) # note that the factor 1/2 is not needed here
    H_z = sum([-m.h_rescaled[i]*Z(i,m.n) for i in range(m.n)])
    H_x = sum([X(i,m.n) for i in range(m.n)])
    
    H = (1-gamma) *(H_zz + H_z) + gamma*H_x
    U = sparse_la.expm(-1j * H * t) # time evoluiton operator
    proposal_mat = np.abs(U)**2

    assert is_stochastic(proposal_mat), 'Proposal matrix is not stochastic.'
    return proposal_mat

def get_proposal_mat_quantum_avg(m, gamma_lims=(0.25, 0.6), gamma_steps=20, time_lims=(2,20)):
    '''Get a quantum proposal matrix for a given Ising model averaged over time and gamma.'''
    # Constructing the Ising and mixing Hamiltonians
    H_zz = sum([-m.J_rescaled[i,j]*Z(i,m.n) @ Z(j,m.n) for i in range(m.n) for j in range(i+1,m.n)]) # note that the factor 1/2 is not needed here
    H_z = sum([-m.h_rescaled[i]*Z(i,m.n) for i in range(m.n)])
    H_x = sum([X(i,m.n) for i in range(m.n)])

    gamma_range = np.linspace(gamma_lims[0], gamma_lims[1], gamma_steps)
    t_min = time_lims[0]
    t_max = time_lims[1]

    d = 2**m.n
    proposal_mat_arr = []
    for gamma in gamma_range:
        H = (1-gamma) *(H_zz + H_z) + gamma*H_x
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
        # Constructing the time-averaged proposal matrix
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

    This function is optimized and is significantly faster than <get_proposal_mat_quantum_avg>.
    ''' 
    def cont_eig(Dlambda):
        t_0, t_f = time_lims[0], time_lims[1] # evolution time t ~ unif(t_0, t_f)
        x = np.sin(Dlambda*t_f) - np.sin(Dlambda*t_0) # from analytical expression for transition probabilities
        return np.divide(2*x/(t_f-t_0), Dlambda, out=np.ones_like(Dlambda), where=(Dlambda!=0) ) 
    
    n = m.n 

    H_zz = sum([-m.J_rescaled[i,j]*Z(i,m.n) @ Z(j,m.n) for i in range(m.n) for j in range(i+1,m.n)]) # note that the factor 1/2 is not needed here
    H_z = H_zz + sum([-m.h_rescaled[i]*Z(i,m.n) for i in range(m.n)])
    H_x = sum([X(i,n) for i in range(n)])

    d = 2**n
    a = np.arange(d**2)
    mask = (a//d >= a%d)
    ones = np.ones(d)

    gamma_range = np.linspace(gamma_lims[0], gamma_lims[1], gamma_steps)

    prop_list = []
    for gamma in gamma_range:
        H = (1-gamma)*H_z + gamma*H_x
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
    return np.all(P >= 0) and np.allclose(col_sums, 1)

def get_transition_matrix(m, T, proposal_mat):
    '''Generate a MC transition probability matrix using Metropolis-Hastings algorithm.'''
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

# Other useful functions
def display_video(video, fps=30):
    """
    Display a 3D numpy array as a video in a Jupyter notebook.

    Parameters:
    array (numpy.ndarray): A 3D numpy array where each slice along the 3rd axis is a frame.
    fps (int): Frames per second, controls the speed of the video.
    """
    n_frames = video.shape[0]
    fig, ax = plt.subplots()
    img = ax.imshow(video[0,:,:], interpolation='nearest', origin='lower')
    for i in range(n_frames):
        frame = video[i,:,:]
        img.set_data(frame)
        display(fig)
        clear_output(wait=True)
        time.sleep(1 / fps)
    plt.close(fig)

