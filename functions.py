import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse.linalg as sparse_la

# Conversion funtions
def int_to_bin(i, n):
    '''Convert a given integer to a bitstring of fixed length using 0 -> 0..00 convention.'''
    return bin(i)[2:].zfill(n)

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

        scale = ( self.n / np.sqrt(la.norm(J, ord='fro')**2 + la.norm(h, ord=2)**2))
        self.J_rescaled = self.J * scale 
        self.h_rescaled = self.h * scale

        self.E = np.zeros((2**self.n))
        for i in range(2**self.n):
            s = int_to_spin(i, self.n)
            self.E[i] = -0.5*(s @ J @ s) - s @ h

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
    J = np.random.standard_normal((n, n))
    J = 0.5 * (J + J.T)  # Ensure symmetry
    np.fill_diagonal(J, 0)  # No self-interactions
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
    p_propose = 1 / sum([scipy.special.binom(m.n,r) for r in range(hamming_radius+1)]) # proposal probability
    proposal_mat = np.zeros((2**m.n, 2**m.n))
    for i in range(2**m.n):
        s_i = int_to_bin(i,m.n)
        for j in range(2**m.n):
            s_j = int_to_bin(j,m.n)
            if hamming(s_i, s_j) > hamming_radius: continue # only local spinflips can be proposed
            proposal_mat[i,j] = p_propose
    return proposal_mat

# Functions for generating a quantum proposal
def get_basis_state(i,n):
    '''Get the computational basis state |i> of n qubits.'''
    assert (i < 2**n) & (i >= 0), 'Bad state.'
    psi = np.zeros(1<<n)
    psi[i] = 1
    return psi

# General helper functions
def is_stochastic(P):
    '''Checks if matrix P is left stochastic, i.e., columns add up to one.'''
    col_sums = np.sum(P, axis=0)
    return np.all(P >= 0) and np.allclose(col_sums, 1)

def get_transition_matrix(m, T, proposal_mat):
    '''Generate a MC transition probability matrix using Metropolis-Hastings algorithm.'''
    P = np.zeros((2**m.n, 2**m.n))
    for i in range(2**m.n):
        for j in range(2**m.n):
            if i == j: continue # sum for j != i
            dE = m.E[i] - m.E[j] # E_new - E_old
            if dE > 0: A = np.exp(-dE / T) # MH acceptance
            else: A = 1 # only compute the exponential if dE > 0 to avoid overflows
            P[i,j] = A * proposal_mat[i,j]
    for i in range(2**m.n):
        P[i,i] = 1 - sum(P[:,i]) # sum for j == i for normalization
    assert is_stochastic(P), 'Something went wrong: P is not stochastic.'
    return P

def get_delta(P):
    '''Calculate the spectral gap of the transition matrix P.'''
    vals = la.eigvals(P) # sparse_la.eigs has problems converging for n > 5
    sorted_vals = sorted(vals, reverse=True)
    assert np.isclose(sorted_vals[0], 1), f'The largest eigenvalue {sorted_vals[0]} differs from 1.'
    return np.real(sorted_vals[0] - sorted_vals[1])