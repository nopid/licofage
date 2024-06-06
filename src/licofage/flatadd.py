from .polytools import companion
from .dfau import dfa
from .intmap import mapper
import numpy as np
from array import array
from . import addcut


def flatten(add, P, cut, verbose):
    n = P.degree()
    M = companion(P).astype(np.intc)
    if verbose:
        print("Expanding DFA...", flush=True, end="")
    add.ssort()
    VF = array("i", [1 if i in add.F else 0 for i in range(add.statecount())])
    VL = np.array([np.array(list(v), dtype=np.intc) for v in add.mv])
    if len(cut) == 1 and n > 1:
        cut = cut * n
    acut = np.array(cut, dtype=np.intc)
    v0 = np.array([0] * n, dtype=np.intc)
    (rn, rq0, RS, RA, RV, RT, RF) = addcut.flatten(
        add.q0, add.S, add.A, add.V, add.T, VF, add.sindex, VL, M, v0, acut, 0
    )
    res = dfa()
    res.q0 = rq0
    res.ms = mapper()
    for i in range(rn):
        res.ms[i]
    res.ma = add.ma.copy()
    res.mv = add.mv.copy()
    res.F = set(RF)
    res.S = RS
    res.A = RA
    res.V = RV
    res.T = RT
    res.ssorted = False
    res.sindex = array("i")
    if verbose:
        print("done.")
    return res
