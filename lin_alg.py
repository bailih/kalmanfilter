import numpy as np

def matrix_square_root(matrix):

    shape = matrix.shape
    if shape[0] != shape[1]:
        raise Exception('Not a square matrix!')

    dim = shape[0]

    eigen_values, P = np.linalg.eig(matrix)

    S = np.zeros((dim, dim))
    for i in range(dim):
        S[i, i] = np.sqrt(eigen_values[i])

    matrix_root = P * S * np.linalg.inv(P)

    return matrix_root

def element_wise_multiplication(matrix1, matrix2):
    shape = matrix1.shape
    if shape!= matrix2.shape:
        raise Exception('Element wise multiplication requires matrices of the same shape.')

    new_matrix = np.mat(np.zeros(shape))

    for i in range(shape[0]):
        for j in range(shape[1]):
            new_matrix[i, j] = matrix1[i, j] * matrix2[i, j]

    return new_matrix

def check_pos_semi_def(matrix, adjust=False):
    '''
        Checks for positive semidefiniteness of a matrix; basically, a symmetric matrix with
        positive eigenvalues.

        If adjust is set to true, this function attempts to force the matrix to be positive
        semidefinite. It does this by replacing any negative eigenvalues with extremely small
        positive ones.
    '''
    shape = matrix.shape
    assert shape[0] is shape[1], 'Not a square matrix, m:{} n:{}'.format(*shape)
    n = shape[0]

    eigen_values, P = np.linalg.eig(matrix)
    try:
        for eig in eigen_values:
            assert eig >= 0, 'Negative eigenvalue found! {}'.format(eig)
        return matrix
    except AssertionError as e:
        if adjust:
            diag = np.mat(np.zeros(shape))
            for i in range(n):
                eig = eigen_values[i]
                diag[i, i] = eig if eig >= 0 else 10**-12
            new_mat = P * diag * np.linalg.inv(P)
            return new_mat
        else:
            raise e


def extend_diagonal(matrix, entry):
    '''
        Takes an n by n matrix and returns an n + 1 by n + 1 matrix with entry as the final
        diagonal.
    '''
    dim = matrix.shape[0]
    new_mat = np.mat(np.zeros((dim + 1, dim + 1)))

    for i in range(dim):
        new_mat[i, i] = matrix[i, i]
    new_mat[dim, dim] = entry

    return new_mat

def extend_vector(vector, entry):
    '''
        Takes an n-dimensional vector and returns a n + 1 dimensional vector with entry as the final
        value.
    '''
    shape = vector.shape
    new_vector = None

    if shape[0] == 1:
        new_vector = np.mat(np.zeros((1, shape[1] + 1)))
        for i in range(shape[1]):
            new_vector[0, i] = vector[0, i]
        new_vector[0, shape[1]] = entry

    elif shape[1] == 1:
        new_vector = np.mat(np.zeros((shape[0] + 1, 1)))
        for i in range(shape[0]):
            new_vector[i, 0] = vector[i, 0]
        new_vector[shape[0], 0] = entry

    else:
        raise Exception('Not a matrix!')

    return new_vector