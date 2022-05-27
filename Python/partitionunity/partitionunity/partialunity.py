import numpy as np
from scipy.spatial import distance_matrix
from scipy.sparse import spdiags
import itertools


# Function for finding the data sites located in each of the q^M blocks
def integer_based_structure(x, q, delta, d):
    """
    Goal: find the data sites located in each of the q^d blocks

    Input   x: nxd numpy array  representing a set of n data sites
            q: number of blocks in one dimension
            delta: radius of PU subdomains
            d: space dimension

    Output: X_block: dictionary {key: values} key represent the index
            of the block  and values is a ndarray that contain indexes
            of points located in k-th block
    """

    X_block = {}
    n, _ = x.shape
    for ind in range(n):
        index = np.ceil(x[ind]/delta)
        index[index == 0] = 1
        index = int(np.sum((index[:-1] - 1)*q**np.arange(d-1, 0, -1)) + index[-1])
        if index in X_block.keys():
            X_block[index] = np.append(X_block[index], [ind], axis=0)
        else:
            X_block[index] = np.array([ind])
    return X_block


# Function that that given a subdomain centre returns the index of
# the square block containing the subdomain centre
def integer_based_containing_query(tilde_x, q, delta, d):
    """
    Goal:   given a subdomain centre returns the index of the
            square block containing the subdomain centre

    Input: tilde_x: subdomain centre
            q: number of blocks in one dimension
            delta: radius of PU subdomains
            d: space dimension

    Output: k: the block containing the subdomain centre
    """
    k = np.ceil(tilde_x/delta)
    k[k == 0] = 1
    k = np.sum((k[:-1] - 1) * q ** np.arange(d - 1, 0, -1)) + k[-1]
    return k


# Function for finding the data sites located in a given subdomain
# and the distances between the subdomain centre and data sites
def integer_based_range_search(tilde_x, delta, x,  X_block, idx_X_neigh_block):
    """
    Goal:   find the data sites located in a given subdomain

    Input:  tilde_x: subdomain centre
            delta: radius of PU subdomain
            x: nxd numpy array  representing a set of n data sites
            idx_X_neigh_block:  dictionary {key: value} key represent
                                the index of the block  and value
                                is a list that contain indexes of
                                points located in k-th block  and in
                                the neighbouring blocks

    Output: n_j:  list of the indexes of the points belonging to a
                  given subdomain
    """

    if idx_X_neigh_block.size != 0:
        n_j = []
        for key in idx_X_neigh_block:
            try:
                for ind in X_block[key]:
                    if np.linalg.norm(tilde_x - x[ind]) <= delta:
                        n_j.append(ind)
            except KeyError:
                pass
        n_j = np.array(n_j)
    else:
        n_j = np.arange(x.shape[0])
    return n_j


def integer_based_neighbourhood(k, q, d):
    """
    Goal: given a block finds neighbouring blocks

    Input:  k: index of the block
            q: number of blocks in one direction
            d: space dimension

    Output: X_NeigBlock: indices of the neighbour blocks
    """

    neigh = np.array([])
    ld = d - 1
    while ld > 0:
        neigh = np.append(neigh, [k + q ** ld, k - q ** ld])
        if ld - 1 > 0:
            neigh = np.append(neigh, np.array([neigh + q ** (ld - 1), neigh - q ** (ld - 1)]))
        ld -= 1

    neigh = np.append(neigh, k)
    neigh = np.append(neigh, [neigh + 1, neigh - 1])
    X_neig_block = neigh[np.logical_and(neigh > 0, neigh <= q ** d)]
    return X_neig_block




def PU(x, f, bar_x, m_d, w,  phi, epsilon):
    """
    Goal: build the partition of unity interpolant and approximate
          the values on a given set of points

    Input:  x: nxd numpy array  representing a set of n data sites
            f: the function values
            bar_x: sxd numpy array  representing a set of s evaluation data
            m_d: number of PU subdomains in one direction
            w: weight function
            phi: radial basis function
            epsilon: the shape parameter

    Output: Pf: sXd numpy array representing the PU fit
    """

    m, d = x.shape
    s, _ = bar_x.shape
    Pf = np.zeros(s)

    # Create m_d^d equally spaced PU centres
    tilde_x = np.array([i for i in itertools.product(np.linspace(0, 1, int(m_d)), repeat=d)])
    # Define the PU radius and the parameter for the weight functions
    radius_par = 2 ** (1 / 2)
    delta = radius_par / m_d
    supp = 1 / delta

    # Parameter for the integer-based partitioning structure
    q = np.ceil(1 / delta)

    # Build the partitioning structure for interpolation and evaluation data
    X_block = integer_based_structure(x, q, delta, d)
    bar_X_block = integer_based_structure(bar_x, q, delta, d)

    # Initialize and compute the Shepard matrix
    sem = w(supp, distance_matrix(bar_x, tilde_x))
    sem = spdiags(1/(sem@np.ones(int(m_d)**d)), 0, s, s)@sem

    # Loop over subdomains
    for j, center in enumerate(tilde_x):

        # Find the box with the j-th PU centre
        k = integer_based_containing_query(center, q, delta, d)
        # Find the interpolation data located in the j-th subdomain
        X_neig_block = integer_based_neighbourhood(k, q, d)
        n_j = integer_based_range_search(center, delta, x, X_block, X_neig_block)

        if n_j.size != 0:
            # Interpolation matrix
            c_j = np.linalg.solve(phi(epsilon, distance_matrix(x[n_j, :], x[n_j, :])), f[n_j])
            bar_X_neig_block = integer_based_neighbourhood(k, q, d)
            s_j = integer_based_range_search(center, delta, bar_x, bar_X_block, bar_X_neig_block)

            if s_j.size != 0:
                # Compute the the local fit
                p_fj = np.dot(phi(epsilon, distance_matrix(bar_x[s_j, :], x[n_j, :])), c_j)
                # Accumulate the global fit
                Pf[s_j] += p_fj * np.array(sem[s_j, j])
    return Pf
