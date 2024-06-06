from .dfau import dfa
from .arraysorter import arraysorter, indexer
from copy import deepcopy


def todfa(q0s, n, S, A, T, F, mapa):
    S = deepcopy(S)
    A = deepcopy(A)
    T = deepcopy(T)
    m = len(S)
    arraysorter([S, A, T])
    sindex = indexer(S, n)
    sindex.append(m - 1)
    seen = set()
    cur = tuple(sorted(list(q0s)))
    seen.add(cur)
    todo = [cur]
    res = dfa()
    res.setinitial(cur)
    while todo:
        cur = todo.pop()
        h = {}
        for s in cur:
            if s in F:
                res.setfinal(cur)
            i = sindex[s]
            if i >= 0:
                while i < m and S[i] == s:
                    a = A[i]
                    ah = h.get(a, set())
                    ah.add(T[i])
                    h[a] = ah
                    i += 1
        for a in h:
            nxt = tuple(sorted(list(h[a])))
            res.addtrans(cur, mapa(a), None, nxt)
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    return res


def rev(aa, follow=None):
    F = [aa.q0]
    if follow:
        todo = [aa.q0]
        seen = set()
        seen.add(aa.q0)
        while todo:
            cur = todo.pop()
            for a, _, t in aa.gettransi(cur):
                if a == follow and t not in seen:
                    todo.add(t)
                    seen.add(t)
                    F.append(t)
    return todfa(aa.F, aa.statecount(), aa.T, aa.A, aa.S, set(F), aa.ma.rev)
