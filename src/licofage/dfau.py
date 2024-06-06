from array import array
from .arraysorter import arraysorter, indexer
from .intmap import mapper, swapper
from . import dfamin
from copy import deepcopy
from itertools import product
import numpy as np


def show(n, m, q0, T, L, H, FF):
    print(f"{n} {m} {q0} {len(FF)}")
    for i in range(m):
        print(f"{T[i]} {L[i]} {H[i]}")
    print(" ".join(str(i) for i in FF))


class dfa:
    def __init__(self):
        self.q0 = None
        self.ms = mapper()
        self.ma = mapper()
        self.mv = mapper()
        self.F = set()
        self.S = array("i")
        self.A = array("i")
        self.V = array("i")
        self.T = array("i")
        self.ssorted = True
        self.sindex = array("i")

    def statecount(self):
        "Count states"
        return len(self.ms)

    def transcount(self):
        "Count transitions"
        return len(self.S)

    def finalcount(self):
        return len(self.F)

    def setinitial(self, s):
        "Set initial state"
        self.q0 = self.ms[s]
        self.ssorted = False

    def getinitial(self):
        return self.ms.rev(self.q0)

    def setfinal(self, s):
        "Set state [s] as final"
        self.F.add(self.ms[s])
        self.ssorted = False

    def isfinal(self, s):
        return self.ms.get(s) in self.F

    def getfinals(self):
        for idxs in self.F:
            yield self.ms.rev(idxs)

    def addtrans(self, s, a, v, t):
        "Add a transition [(s,a,v,t)] -- NO COLLISION CHECK"
        self.S.append(self.ms[s])
        self.A.append(self.ma[a])
        self.V.append(self.mv[v])
        self.T.append(self.ms[t])
        self.ssorted = False

    def ssort(self):
        "Sort transitions by initial state"
        if self.ssorted:
            return
        arraysorter([self.S, self.A, self.V, self.T])
        self.sindex = indexer(self.S, self.statecount())
        self.ssorted = True

    def copy(self):
        res = dfa()
        res.q0 = self.q0
        res.ms = self.ms.copy()
        res.ma = self.ma.copy()
        res.mv = self.mv.copy()
        res.F = self.F.copy()
        res.S = deepcopy(self.S)
        res.A = deepcopy(self.A)
        res.V = deepcopy(self.V)
        res.T = deepcopy(self.T)
        res.ssorted = self.ssorted
        res.sindex = deepcopy(self.sindex)
        return res

    def alltransi(self):
        "Generate indexes of transition (s,a,v,t)"
        for i in range(len(self.S)):
            yield (self.S[i], self.A[i], self.V[i], self.T[i])

    def alltrans(self, sorted=False):
        "Generate all transitions"
        if sorted:
            self.ssort()
        for s, a, v, t in self.alltransi():
            yield (self.ms.rev(s), self.ma.rev(a), self.mv.rev(v), self.ms.rev(t))

    def gettransi(self, idxs):
        "Generate all transitions from state of index [idxs]"
        self.ssort()
        m = self.transcount()
        i = self.sindex[idxs]
        if i >= 0:
            while i < m and self.S[i] == idxs:
                yield (self.A[i], self.V[i], self.T[i])
                i += 1

    def gettrans(self, s):
        "Generate all transitions from states [s]"
        idxs = self.ms.get(s)
        for a, v, t in self.gettransi(idxs):
            yield (self.ma.rev(a), self.mv.rev(v), self.ms.rev(t))

    def delta(self, s, v0, red, *l):
        self.ssort()
        m = self.transcount()
        idxs = self.ms.get(s)
        cur = idxs
        vcur = v0
        for a in l:
            idxa = self.ma.get(a)
            i = self.sindex[cur]
            if i < 0:
                return None
            while i < m and self.S[i] == idxs and self.A[i] != idxa:
                i += 1
            if i == m or self.A[i] != idxa:
                return None
            v = self.mv.rev(self.V[i])
            cur = self.T[i]
            vcur = red(vcur, v)
        return vcur, cur

    def maptrans(self, f):
        "Remap [a] and [v] from every transition [(s,a,v,t)] using [(b,w)=f(a,v)] -- NO COLLISION CHECK"
        res = dfa()
        res.q0 = self.q0
        res.ms = self.ms.copy()
        res.ma = mapper()
        res.mv = mapper()
        res.F = self.F.copy()
        res.S = array("i")
        res.A = array("i")
        res.V = array("i")
        res.T = array("i")
        res.ssorted = False
        res.sindex = array("i")
        for s, a, v, t in self.alltrans():
            b, w = f(a, v)
            res.addtrans(s, b, w, t)
        return res

    def trim(self):
        "Trim automaton"
        res = dfa()
        seen = set()
        res.setinitial(self.ms.rev(self.q0))
        back = {}
        todo = [self.q0]
        seen.add(self.q0)
        while todo:
            si = todo.pop()
            for ai, vi, ti in self.gettransi(si):
                lti = back.get(ti, [])
                lti.append((si, ai, vi))
                back[ti] = lti
                if ti not in seen:
                    seen.add(ti)
                    todo.append(ti)
        coacc = set()
        todo = []
        for si in self.F:
            if si in back:
                s = self.ms.rev(si)
                res.setfinal(s)
                coacc.add(si)
                todo.append(si)
        while todo:
            ti = todo.pop()
            t = self.ms.rev(ti)
            for si, ai, vi in back[ti]:
                s = self.ms.rev(si)
                a = self.ma.rev(ai)
                v = self.mv.rev(vi)
                res.addtrans(s, a, v, t)
                if si not in coacc:
                    coacc.add(si)
                    todo.append(si)
        return res

    def minimize(self):
        "Minimize [a] loosing the vectors (copy of labels from now)"
        n = self.statecount()
        m = self.transcount()
        res = dfa()
        (n, m, q0, T, L, H, FF) = dfamin.minimize(
            n,
            m,
            self.q0,
            deepcopy(self.S),
            deepcopy(self.A),
            deepcopy(self.T),
            array("i", sorted(list(self.F))),
        )
        res.q0 = q0
        res.ms = mapper()
        for i in range(n):
            res.ms[i]
        res.ma = self.ma.copy()
        res.mv = self.ma.copy()
        res.F = set(FF)
        res.S = T
        res.A = L
        res.V = deepcopy(L)
        res.T = H
        res.ssorted = False
        res.sindex = array("i")
        return res

    def todot(self, f, vformat=str, vsep=":"):
        """Write dot representation into [f] object using [vformat] to
        format vector component. If [vformat] is [None] or return
        [None], the vector component is ignored, otherwise it
        appears on the edge separated from the label by [vsep]"""
        f.write("""digraph G {\nrankdir=LR\nnode [shape=circle]\n"" [shape=none]\n""")
        for sidx in self.F:
            f.write(f'"{self.ms.rev(sidx)}" [shape=doublecircle]\n')
        f.write(f'"" -> "{self.ms.rev(self.q0)}"\n')
        h = {}
        for s, a, v, t in self.alltrans():
            r = vformat(v) if vformat is not None else None
            if r is None:
                lbl = str(a)
            else:
                lbl = f"{a}{vsep}{r}"
            ll = h.get((s, t), [])
            ll.append(lbl)
            h[(s, t)] = ll
        nl = "\n"
        for (s, t), lbl in h.items():
            f.write(f'"{s}" -> "{t}" [label="{nl.join(lbl)}"]\n')
        f.write("}\n")

    def inferalpha(self):
        "Attempt to reconstruct DFA alphabet from transitions."
        v = next(iter(self.ma))
        if type(v) is not tuple:
            return [sorted(list(self.ma))]
        n = len(v)
        s = [set() for _ in range(n)]
        for x in self.ma:
            for i, v in enumerate(x):
                s[i].add(v)
        return [sorted(list(x)) for x in s]

    def toWalnut(self, f, alpha=None, outmap=None):
        """Write Walnut representation into [f] object using
        only [(s,a,t)] labels. If [outmap] is given, use
        it to determine output of DFAO. Otherwise use
        1 for accepting states. Alphabets are given by
        [alpha] string if given. Otherwise alphabets are
        inferred from labels."""
        self.ssort()
        swp = swapper(0, self.q0)
        if alpha is None:
            alpha = " ".join(
                "{" + ", ".join(map(str, x)) + "}" for x in self.inferalpha()
            )
        f.write(f"{alpha}\n")
        for idxs in range(self.statecount()):
            ridxs = swp[idxs]
            if outmap is None:
                v = 1 if ridxs in self.F else 0
            else:
                v = outmap(self.ms.rev(ridxs))
            f.write(f"\n{idxs} {v}\n")
            for idxa, _, ridxt in self.gettransi(ridxs):
                a = self.ma.rev(idxa)
                idxt = swp[ridxt]
                f.write(
                    f"{' '.join(map(str,a)) if type(a) is tuple else str(a)} -> {idxt}\n"
                )


def prod(f, *aa):
    res = dfa()
    res.setinitial(tuple(map(lambda a: a.getinitial(), aa)))
    for x in set(product(*map(lambda a: a.getfinals(), aa))):
        res.setfinal(x)
    for tl in product(*map(lambda a: a.alltrans(), aa)):
        sl = tuple(map(lambda x: x[0], tl))
        al = tuple(map(lambda x: x[1], tl))
        vl = tuple(map(lambda x: x[2], tl))
        tl = tuple(map(lambda x: x[3], tl))
        b, w = f(*al, *vl)
        res.addtrans(sl, b, w, tl)
    return res


def lincomb(a, coef, objv):
    "linear combinations of automaton [a] : \sum_i coef[i] * x[i] = objv"
    nc = len(coef)

    def aux(*p):
        nonlocal nc, coef
        ra = tuple(p[:nc])
        rv = tuple(sum(ci * np.array(xi) for ci, xi in zip(coef, p[nc:])))
        return (ra, rv)

    return prod(aux, *([a] * len(coef)))


def addition(a):
    return lincomb(a, [1, 1, -1], 0)
