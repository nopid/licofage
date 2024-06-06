import numpy as np
from .polytools import companion
from . import matrix


class rebaser:
    def __init__(self, P, u):
        n = P.degree()
        M = companion(P)
        VV = np.stack([np.dot(np.linalg.matrix_power(M, i), u) for i in range(n)])
        self.mat = matrix.inv(VV)
        assert self.mat is not None, f"bad vector {u} with {P} :\n{VV}"

    def topoly(self, v):
        return np.dot(self.mat, v)
