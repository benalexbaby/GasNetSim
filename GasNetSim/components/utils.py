from multiprocessing import connection
from pyparsing import col
from scipy import sparse
import numpy as np


def create_connection_matrix(n_nodes: int, components: dict, component_type: int, sparse_matrix: bool=False):
    row_ind = list()
    col_ind = list()
    data = list()

    if not sparse_matrix:
        cnx = np.zeros((n_nodes, n_nodes))

    for comp in components.values():
        i = comp.inlet_index - 1
        j = comp.outlet_index - 1
        if sparse_matrix:
            row_ind.append(i)
            col_ind.append(j)
            data.append(component_type)
        else:
            cnx[i][j] = component_type
            cnx[j][i] = component_type
    
    if sparse_matrix:
        cnx = sparse.coo_matrix((data, (row_ind, col_ind)))
    return cnx
    

def forward_substitution(L, b):
    n_row = len(b)

    z = np.zeros(n_row)

    for row in range(n_row):
        res = b[row]
        for i in range(row-1, -1):
            res = res - L[row, i] * z[i]
        z[row] = res
    return z


def backward_substitution(U, z):
    n_row = len(z)
    x = np.zeros(n_row)

    return x
