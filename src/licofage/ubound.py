import numpy as np
from math import floor


def row(i, r, n):
    if i < len(r):
        return [pow(r[i], j) for j in range(n)]
    else:
        return [1 if j == i - len(r) else 0 for j in range(n)]


def matP(P):
    n = P.degree()
    r = list(P.roots())
    r.sort(key=lambda v: abs(v), reverse=True)
    assert len(r) == n
    k = 0
    while abs(r[-1]) < 1e-4:
        k += 1
        r.pop(-1)
    return np.array([row(i, r, n) for i in range(n)]).transpose(), abs(r[0]), k


def boundit(P, u, vobj):
    k = P.degree()
    A, th, kX = matP(P)
    gammas = np.dot(np.linalg.inv(A), u[:k])
    alpha = sum([abs(gammas[i]) / (1 - abs(A[1, i])) for i in range(1, k - kX)])
    alpha += sum([abs(gammas[i]) for i in range(k - kX, k)])
    beta = alpha + abs(gammas[0]) * th / (th - 1)
    K = abs(vobj)
    return [(K * pow(th, i)) for i in range(k)], [
        (alpha + beta * pow(th, i)) for i in range(k)
    ]
    return [floor(K * pow(th, i)) for i in range(k)], [
        floor(alpha + beta * pow(th, i)) for i in range(k)
    ]
