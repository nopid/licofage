"""Precise (but slow) operations on integer matrices"""
import numpy as np
from fractions import Fraction
from .bareiss import adjugate


def inv(M):
    s, d, A = adjugate(M)
    if A is None:
        return None
    R = Fraction(s, d) * np.array(A)
    return R
