import numpy as np
import scipy as sc
from scipy.sparse import coo_matrix, hstack, csr_matrix, bmat, dia
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import eigs
import time


class Laplacian:
    """
        The Laplacian class implements the Laplace Operator for the given problem.

        For now, I only implemented the Neuman BC with this class

        Attributes:
            params (tuple): a tuple containing the discretization steps of the grid
        """

    def __init__(self, params):
        self.NX, self.NY = params
        self._construct_sparse_neuman()

    def _construct_sparse_neuman(self):
        #for square grid this is the same
        L = sc.sparse.eye(self.NX, self.NX- 1)
        L.setdiag(-1, k = -1)
        I = sc.sparse.eye(self.NX, self.NX)
        Ax = L.dot(L.T)
        self.A = sc.sparse.kron(I, Ax) + sc.sparse.kron(Ax, I)
