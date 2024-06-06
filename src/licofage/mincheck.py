from .polytools import companion
import numpy as np


def revector(fa, lcfa, P):
    "try to associate a single vector to every state of fa"
    ok = True
    n = P.degree()
    M = companion(P)
    vecmap = {}
    mq0 = fa.getinitial()
    lq0 = lcfa.getinitial()
    v0 = tuple([0] * n)
    z = (mq0, lq0, v0)
    seen = set()
    todo = [z]
    seen.add(z)
    while todo:
        (ms, ls, vs) = todo.pop()
        if ms in vecmap and vs not in vecmap[ms]:
            ok = False
        daset = vecmap.get(ms, set())
        daset.add(vs)
        vecmap[ms] = daset
        tl = {}
        tm = {}
        for a, v, lt in lcfa.gettrans(ls):
            tl[a] = (tuple(np.dot(M, vs) + v), lt)
        for a, _, mt in fa.gettrans(ms):
            tm[a] = mt
        for a in tm:
            (vt, lt) = tl[a]
            mt = tm[a]
            z = (mt, lt, vt)
            if z not in seen:
                seen.add(z)
                todo.append(z)
    if ok:
        print(f"Re-vectorization ok! ({len(vecmap)} states)")
    else:
        print("#### Re-vectorization KO")
        for q, s in vecmap.items():
            if len(s) > 1:
                print(f"{q:6d} : {s}")
