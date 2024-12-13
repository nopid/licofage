"""A convenient lightweight library to construct and combine sequence automata."""

import numpy as np
from .misc import poly2str, str2poly
from .subst import (
    str2sub,
    sub2str,
    addpoly,
    substdfa,
    substparentdfa,
    dfasubst,
    substpoly,
    alpha,
    alphaim,
)
from .subst import block as doblock
from .polytools import (
    gen_pisot_root,
    companion,
    x as poly_x,
    one as poly_on,
    one as poly_one,
    naivelcm,
)
from .polybound import compcut, realcut
from .flatadd import flatten
from .dfau import prod, dfa
from pathlib import Path as P

verbose = False
stats = False
destdir = "/tmp/out"


def setparams(v, s, d):
    global verbose, stats, destdir
    verbose = v
    stats = s
    destdir = P(d)


def writefile(name, data, dir=P("Command Files")):
    genparentdir(destdir / dir / P(f"{name}.txt")).write_text(data)


def genparentdir(fname):
    parent = fname.parent
    if not parent.is_dir():
        parent.mkdir(parents=True)
    return fname


def prstats(title, fa):
    print(
        f"{title}: {fa.statecount()} states, {fa.finalcount()} final states, {fa.transcount()} transitions."
    )


def unpoly(e):
    if all(x >= 0 for x in e):
        return parry(*e)
    n = len(e)
    d = [0] * n
    for m in range(n):
        for i in range(m - 1):
            d[i] = e[i]
        d[m - 1] = e[m - 1] - 1
        for i in range(n - m):
            d[m + i] = d[i] + e[m + i]
        if all(x >= 0 for x in d):
            return parry(*d[: n - m]).tail(*d[n - m :])
    return None


class poly:
    def __init__(self, s):
        self.P = str2poly(s)

    def __repr__(self):
        return f"poly({poly2str(self.P)})"

    def image(self, s):
        return morphic(self, s)

    def subst(self):
        return unpoly([-int(x) for x in list(self.P)[::-1][1:]]).subst()


class parry:
    def __init__(self, *coeffs):
        self.coeffs = coeffs

    def __repr__(self):
        return f"parry({str(self.coeffs)[1:-1]})"

    def tail(self, *coeffs):
        return parrytail(self.coeffs, coeffs)

    def image(self, s):
        return morphic(self, s)

    def subst(self):
        h = dict()
        m = len(self.coeffs) - 1
        for i, x in enumerate(self.coeffs):
            h[str(i)] = "0" * x + (str(i + 1) if i < m else "")
        return ("0", h)


class parrytail:
    def __init__(self, u, v):
        self.u = u
        self.v = v

    def __repr__(self):
        return f"parry({str(self.u)[1:-1]}[{str(self.v)[1:-1]}]*)"

    def image(self, s):
        return morphic(self, s)

    def subst(self):
        h = dict()
        coeffs = self.u + self.v
        k = len(self.u)
        m = len(coeffs) - 1
        for i, x in enumerate(coeffs):
            h[str(i)] = "0" * x + (str(i + 1) if i < m else str(k))
        return ("0", h)


class subst:
    def __init__(self, s):
        self.a = s.strip()[0]
        self.h = str2sub(s)

    def __repr__(self):
        return f"subst({sub2str(self.h)})"

    def image(self, s):
        return morphic(self, s)

    def subst(self):
        return (self.a, self.h)


class morphic:
    def __init__(self, base, s):
        self.base = base
        self.h = str2sub(s)

    def __repr__(self):
        return f"morphic({repr(self.base)} | {sub2str(self.h)})"

    def image(self, s):
        f = str2sub(s)
        h = dict()
        for x, u in self.h.items():
            h[x] = "".join([f[y] for y in u])
        return morphic(self.base, sub2str(h))

    def subst(self):
        return self.base.subst()

    def morphism(self):
        return self.h


class block:
    def __init__(self, base, n=2):
        self.base = base
        self.n = n

    def __repr__(self):
        return f"block({repr(self.base)}, {self.n})"

    def image(self, s):
        return morphic(self, s)

    def subst(self):
        return doblock(*self.base.subst(), self.n)


class address:
    def __init__(self, base, nsname, **weights):
        self.base = base
        self.ns = nsname
        self.w = weights
        if verbose:
            print(f">>> Address automaton for {repr(self.base)}")
        self.a, self.fh = self.base.subst()
        alph = alpha(self.fh)
        morphic = hasattr(self.base, "morphism")
        self.gh = self.base.morphism() if morphic else None
        self.alpha = alphaim(self.gh) if morphic else alph
        self.P = poly_one
        iw = None
        self.rw = {x: 1 for x in self.alpha}
        if morphic:
            iw = [len(self.gh[x]) for x in alph]
            if max(iw) > 1:
                self.P *= poly_x
        if weights:
            rw = dict()
            dw = weights.get("_", 1)
            for x in self.alpha:
                rw[x] = weights.get(x, dw)
            self.rw = rw
            if morphic:
                iw = [sum([rw[y] for y in self.gh[x]]) for x in alph]
            else:
                iw = [rw[x] for x in alph]
        self.iw = iw
        self.P *= addpoly(self.fh, iw)
        self.n = self.P.degree()
        if verbose:
            print(f"  DT polynomial: {poly2str(self.P)} (n={self.n})")
        self.r = gen_pisot_root(self.P)
        if self.r is None and verbose:
            print("  *** WARNING *** the substitution is not (X^k-)Pisot!")
        elif self.r is not None and verbose:
            kX, th = self.r
            print(f"  θ={th}" + (f" (X^{kX})" if kX > 0 else ""))
        self.fa = self.gen_fa(self.n)
        if stats:
            prstats("  DT automaton", self.fa)
        if verbose:
            print()

    def gen_fa(self, n):
        if n < self.n:
            raise ValueError("Shorter vectors than base!")
        fa = substdfa(self.fh, self.a, n, self.iw)
        if self.gh is not None:
            newfa = dfa()
            alen = [0] * len(fa.ma)
            for s, a, v, t in fa.alltrans(True):
                k = len(self.gh[t])
                alen[a] = max(alen[a], k)
            aidx = []
            cur = 0
            for x in alen:
                aidx.append(cur)
                cur += x
            for i in range(sum(alen)):
                newfa.ma[i]
            newfa.setinitial(fa.getinitial())
            zv = [0] * (self.n - 1)
            for s, a, v, t in fa.alltrans(True):
                u = self.gh[t]
                k = len(u)
                if k > 0:
                    newfa.setfinal(t)
                ai = aidx[a]
                newfa.addtrans(s, ai, v, t)
                iv = v[0]
                for i, y in list(enumerate(u))[1:]:
                    di = sum(self.rw[y] for y in u[:i])
                    newfa.addtrans(s, ai + i, tuple([iv + di] + zv), y)
                    newfa.setfinal(y)
            fa = newfa
        return fa

    def state_output(self, x):
        if self.gh is None:
            return x
        if x in self.fh:
            return self.gh[x][0]
        return x

    def __repr__(self):
        return f"address({repr(self.base)}, '{self.ns}'{', ' + str(self.w)[1:-1] if self.w else ''})"

    def poly(self):
        return self.P

    def gen_ns(self, format="Walnut", minimize=True, cut=None):
        if verbose:
            print(f">>> {repr(self)} is generating it's NS")
        fa = self.fa
        if minimize:
            fa = fa.minimize()
        if stats:
            prstats("  Numeration system automaton", fa)
        fname = genparentdir(destdir / P("Custom Bases") / P(f"msd_{self.ns}.txt"))
        with open(fname, "w") as f:
            print(f"  Writing {fname} in format {format}")
            if format == "Walnut":
                fa.toWalnut(f)
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")
        (self + self - self).gen_dfa(
            f"msd_{self.ns}_addition",
            cut=cut,
            dir=P("Custom Bases"),
            set_ns=False,
            preverb="  ",
        )
        if verbose:
            print()

    def gen_word_automaton(self, name=None, format="Walnut"):
        if name is None:
            name = self.ns.title()
        if verbose:
            print(
                f">>> {repr(self)} is generating it's word automaton under name {name}"
            )
        fa = self.fa
        fname = genparentdir(
            destdir / P("Word Automata Library") / P(f"{name.title()}.txt")
        )
        with open(fname, "w") as f:
            if verbose:
                print(f"  Writing {fname} in format {format}\n")
            if format == "Walnut":
                fa.toWalnut(
                    f,
                    alpha=f"msd_{self.ns}",
                    outmap=lambda x: self.alpha.index(self.state_output(x)),
                )
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")
        if verbose:
            print()

    def __mul__(self, other):
        assert isinstance(other, int)
        return lincomb((other, self))

    def __rmul__(self, other):
        assert isinstance(other, int)
        return lincomb((other, self))

    def __add__(self, other):
        return lincomb((1, self)) + other

    def __radd__(self, other):
        return other + lincomb((1, self))

    def __sub__(self, other):
        return lincomb((1, self)) - other

    def __rsub__(self, other):
        return other - lincomb((1, self))


def opposite(xs):
    return [(-v, y) for (v, y) in xs]


class lincomb:
    def __init__(self, *xs):
        self.xs = xs

    def __add__(self, other):
        if type(other) == lincomb:
            return lincomb(*self.xs, *other.xs)
        else:
            return lincomb(*self.xs, (1, other))

    def __radd__(self, other):
        if type(other) == lincomb:
            return lincomb(*other.xs, *self.xs)
        else:
            return lincomb((1, other), *self.xs)

    def __sub__(self, other):
        if type(other) == lincomb:
            return lincomb(*self.xs, *opposite(other.xs))
        else:
            return lincomb(*self.xs, (-1, other))

    def __rsub__(self, other):
        if type(other) == lincomb:
            return lincomb(*other.xs, *opposite(self.xs))
        else:
            return lincomb((1, other), *opposite(self.xs))

    def __repr__(self):
        return "".join(
            map(
                lambda z: f"{'+' if z[0]>=0 else '-'}{(str(abs(z[0]))+'*') if abs(z[0])!=1 else ''}{repr(z[1])}",
                self.xs,
            )
        )

    def compute_fa(self):
        self.compute_poly()
        nc = len(self.xs)
        coef, addas = zip(*self.xs)
        fas = map(lambda x: x.gen_fa(self.n), addas)
        alpha = " ".join(map(lambda x: f"msd_{x.ns}", addas))

        def aux(*p):
            nonlocal nc, coef
            ra = tuple(p[:nc])
            rv = tuple(sum(ci * np.array(xi) for ci, xi in zip(coef, p[nc:])))
            return (ra, rv)

        return (prod(aux, *fas), alpha)

    def compute_poly(self):
        self.P = poly_one
        _, addas = zip(*self.xs)
        for adda in addas:
            self.P = naivelcm(self.P, adda.poly())
        self.n = self.P.degree()
        if verbose:
            print(f"  Combination polynomial is {poly2str(self.P)}")
        self.r = gen_pisot_root(self.P)
        if self.r is None and verbose:
            print("  *** WARNING *** the substitution is not (X^k-)Pisot!")
        elif self.r is not None and verbose:
            kX, th = self.r
            print(f"  θ={th}" + (f" (X^{kX})" if kX > 0 else ""))

    def gen_dfa(
        self,
        name,
        cut=None,
        vobj=0,
        vector=None,
        minimize=True,
        dir=P("Automata Library"),
        set_ns=True,
        format="Walnut",
        preverb=">>> ",
        debug=None,
    ):
        if verbose:
            print(f"{preverb}{repr(self)} is generating it's dfa under name {name}")
        fname = genparentdir(destdir / dir / P(f"{name}.txt"))
        lcfa, alpha = self.compute_fa()
        if debug:
            with open(genparentdir(destdir / P("debug") / P(f"{name}.dot")), "w") as f:
                lcfa.todot(f)
        if vector is None:
            vector = [0] * (self.n - 1) + [1]
        if stats:
            prstats("  Combination automaton", lcfa)
        if cut is None:
            if self.r is None:
                raise ValueError("the substitution is not (X^k-)Pisot!")
            cut = compcut(lcfa, self.P, vector, 0, verbose, preverb="  ")
        elif verbose:
            print(f"  Bypassing bound computation! cut={cut}")
        fa = flatten(lcfa, self.P, cut, verbose)
        if stats:
            prstats(f"  Flattened automaton ", fa)
        fa = fa.trim()
        if stats:
            prstats(f"  Trimmed flattened automaton ", fa)
        if verbose:
            obscut = realcut(fa, companion(self.P))
            print(f"  Real bound is {obscut} (vs {cut}).")
        if minimize:
            fa = fa.minimize()
            if stats:
                prstats("  After minimization", fa)
        with open(fname, "w") as f:
            if verbose:
                print(f"  Writing {fname} in format {format}\n")
            if format == "Walnut":
                fa.toWalnut(f, alpha if set_ns else None)
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")
