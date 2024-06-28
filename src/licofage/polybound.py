from .polytools import unfrac
from .lrseq import rebaser
from .ubound import boundit
from .misc import stra
import numpy as np


def clean(l, vl, s):
    return (l, tuple(map(lambda v: tuple(list(v)), vl)), s)


def value(vl):
    res = 0
    for i, l in enumerate(vl):
        res += l[i]
    return res


def compcut(add, P, U, vobj, verbose, preverb=""):
    n = P.degree()
    boundK, boundC = boundit(unfrac(P), U, vobj)
    if verbose:
        print(f"{preverb}(P,u) vector bound: {stra(boundK)} + {stra(boundC)}*C")
    reb = rebaser(P, U)
    addpo = add.maptrans(lambda a, v: (a, tuple(reb.topoly(v))))
    mi, ma = [0] * n, [0] * n
    for s, a, v, t in addpo.alltrans():
        for i in range(n):
            if v[i] < mi[i]:
                mi[i] = v[i]
            if v[i] > ma[i]:
                ma[i] = v[i]
    C = max(abs(sum(mi)), abs(sum(ma)))
    cut = [b + C * a for a, b in zip(boundC, boundK)]
    if verbose:
        print(f"{preverb}Fast poly bound: C in [{stra(mi)}, {stra(ma)}] => {stra(cut)}")
    return cut


def realcut(res, M):
    n = len(next(res.alltrans())[2])
    obscut = [0] * n
    q0 = res.getinitial()
    seen = set()
    todo = [(q0, tuple([0] * n))]
    while todo:
        cur, curv = todo.pop()
        for i in range(n):
            obscut[i] = max(obscut[i], abs(curv[i]))
        for a, v, t in res.gettrans(cur):
            if t not in seen:
                seen.add(t)
                nxtv = tuple(list(np.dot(M, curv) + np.array(v)))
                todo.append((t, nxtv))
    return obscut
