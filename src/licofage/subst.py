import numpy as np
from .polytools import commonpoly
from .dfau import dfa
from .intmap import mapper

ALPHABET = "abcdefghijklmnopqrstuvwxyz"


def str2sub(s):
    res = {}
    if "->" in s:
        for l in s.strip().split(","):
            (a, u) = l.strip().split("->")
            res[a.strip()] = u.strip()
    else:
        for i, u in enumerate(s.strip().split("/")):
            res[str(i)] = u.strip()
    return res


def sub2str(h):
    return ", ".join(f"{a}->{sa}" for a, sa in sorted(h.items()))


def canon(s, a, alpha=ALPHABET):
    h = str2sub(s)
    res = []
    ms = mapper()
    seen = set()
    todo = [a]
    seen.add(a)
    while todo:
        cur = todo.pop()
        s = f"{ alpha[ms[cur]] }->"
        for x in h[cur]:
            s += alpha[ms[x]]
            if x not in seen:
                seen.add(x)
                todo.append(x)
        res.append(s)
    return ", ".join(sorted(res))


def alpha(s):
    return sorted(list(s.keys()))


def alphaim(s):
    return sorted(list(set(sum(map(list, s.values()), []))))


def seqalpha(s):
    k = 0
    for _, v in s.items():
        k = max(k, len(v))
    return list(range(k))


def substmat(s):
    "compute the adjacency matrix of substitution [s]"
    alph = alpha(s)
    return np.array([[s[i].count(j) for j in alph] for i in alph])


def substpoly(s):
    return np.polynomial.Polynomial(np.poly(substmat(s))[::-1])


def seq(s, iw=None):
    """compute the [2k] first terms of the length sequences for every letters of substitution [s]
    starting from initial weight values [iw]"""
    alph = alpha(s)
    A = substmat(s)
    k = 2 * len(alph)
    vv = []
    if iw is None:
        curv = np.array([[1]] * len(alph))
    else:
        curv = np.array([[i] for i in iw])
    for _ in range(k):
        vv.append(curv)
        curv = np.dot(A, curv)
    res = {x: [] for x in alph}
    res[""] = [0] * k
    for v in vv:
        for x, k in zip(alph, v):
            res[x].append(k[0])
    for x in res:
        res[x] = np.array(res[x])
    return res


def combiseq(bs, u):
    "compute the first terms of word [u] for a dictionnary of sequences [bs]"
    if u == "":
        return bs[""]
    return list(sum(bs[x] for x in u))


def addpoly(s, iw=None):
    "compute the minimal reccurence polynomial for addition of [s]"
    vv = seq(s, iw)
    pp = []
    for x in s:
        u = s[x]
        for i in range(len(u) - 1):
            v = u[: i + 1]
            pp.append(combiseq(vv, v))
    return commonpoly(*pp)


def substdfa(s, a, k=0, iw=None):
    """compute the sequence value automaton for [s] from [a] with vectors of size [k]
    starting with inital weight values [iw]"""
    assert (
        s[a][0] == a
    ), f'the image of "{a}" does not starts with "{a}" in the substitution'
    res = dfa()
    vv = seq(s, iw)
    for x in s:
        u = s[x]
        for i in range(len(u)):
            v = combiseq(vv, u[:i])[:k]
            assert (
                len(v) == k
            ), f"oups... something strange with sequences lengths ({k}, {v})"
            y = u[i]
            res.addtrans(x, i, tuple(v), y)
    res.setinitial(a)
    for x in s:
        res.setfinal(x)
    return res


def substparentdfa(s, a, k=0):
    "compute the parent sequence value automaton for [s] from [a] with vectors of size [k]"
    cla = substdfa(s, a, k)
    res = dfa()
    oq0 = cla.getinitial()
    q0 = (oq0, (oq0, 0))
    res.setinitial(q0)
    seen = set()
    todo = [q0]
    seen.add(q0)
    while todo:
        ns = todo.pop()
        (s, orig) = ns
        if cla.isfinal(s):
            res.setfinal(ns)
        for a, v, t in cla.gettrans(s):
            nt = (t, (s, a))
            res.addtrans(ns, a, v, nt)
            if nt not in seen:
                seen.add(nt)
                todo.append(nt)
    return res


def dfasubst(aa):
    "compute the substitution of a prefix-closed automaton [aa]"
    h = {}
    alpha = ALPHABET
    ms = mapper()
    ms[aa.getinitial()]
    for ss, a, _, tt in aa.alltrans():
        s = alpha[ms[ss]]
        t = alpha[ms[tt]]
        if not aa.isfinal(ss):
            return None
        if type(a) is not int:
            return None
        ls = h.get(s, [])
        ls.append((a, t))
        h[s] = ls
    for s in h:
        ls = sorted(h[s])
        fs = ""
        for i, (j, c) in enumerate(ls):
            if i != j:
                return None
            fs += c
        h[s] = fs
    return h
