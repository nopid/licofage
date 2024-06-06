from .subst import (
    str2sub,
    sub2str,
    addpoly,
    substdfa,
    substparentdfa,
    dfasubst,
    substpoly,
    alpha,
)
from .polytools import gen_pisot_root, companion
from .polybound import compcut, realcut
from .flatadd import flatten
from .dfau import lincomb
from .misc import poly2str
from .intmap import mapper
from .lrseq import rebaser
from .mincheck import revector
from . import nfa
from numpy.polynomial import Polynomial as Poly


def prstats(title, fa):
    print(
        f"{title}: {fa.statecount()} states, {fa.finalcount()} final states, {fa.transcount()} transitions."
    )


def img(a, h, n):
    cur = a
    while len(cur) < n:
        nxt = ""
        for c in cur:
            nxt += h[c]
        cur = nxt
    return cur[:n]


class Study:
    def __init__(self, s, verbose=False, stats=False, simplify=True):
        "An API to study a given textual substitution s"
        self.verbose = verbose
        self.stats = stats
        if verbose:
            print(f"Studying substitution {s}")
        self.a = s.strip()[0]
        self.h = str2sub(s)
        self.origh = self.h
        self.origa = self.a
        if verbose:
            print(img(self.a, self.h, 70) + "...")
            print(f"Substitution polynomial: {poly2str(substpoly(self.h))}")
        if simplify:
            self.h = dfasubst(substdfa(self.h, self.a).minimize())
            self.a = "a"
            if verbose:
                print(f"Minimized substitution {sub2str(self.h)}")
                print(f"Substitution polynomial: {poly2str(substpoly(self.h))}")
        self.P = addpoly(self.h)
        self.n = self.P.degree()
        if verbose:
            print(f"Addition polynomial: {poly2str(self.P)} (n={self.n})")
        self.r = gen_pisot_root(self.P)
        if self.r is None and verbose:
            print("*** WARNING *** the substitution is not (X^k-)Pisot!")
        elif self.r is not None and verbose:
            kX, th = self.r
            print(f"θ={th}" + (f" (X^{kX})" if kX > 0 else ""))
        self.fa = substdfa(self.h, self.a, self.n)
        if self.stats:
            prstats("Substitution automaton", self.fa)

    def gendfao(self, fname, format):
        fa = substdfa(self.origh, self.origa)
        with open(fname, "w") as f:
            print(f">>> Writing {fname} in format {format}\n")
            if format == "Walnut":
                fa.toWalnut(f, outmap=lambda x: alpha(self.origh).index(x))
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")

    def genparentdfao(self, fname, format):
        fa = substparentdfa(self.origh, self.origa)
        out = mapper()
        with open(fname, "w") as f:
            print(f">>> Writing {fname} in format {format}\n")
            if format == "Walnut":
                fa.toWalnut(f, outmap=lambda x: out[x[1]])
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")

    def gennumsys(self, fname, msd, format, minimize):
        fa = self.fa
        if not msd:
            fa = nfa.rev(fa)
        if minimize:
            fa = fa.minimize()
        if self.stats:
            prstats("Numeration system automaton", fa)
        with open(fname, "w") as f:
            print(f">>> Writing {fname} in format {format}\n")
            if format == "Walnut":
                fa.toWalnut(f)
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")

    def genlinear(
        self,
        fname,
        lineq,
        msd,
        format,
        minimize,
        seq,
        poly,
        vector,
        cut,
        dorevector=False,
    ):
        if vector is None:
            vector = [0] * (self.n - 1) + [1]
        lcfa = lincomb(self.fa, *lineq)
        if self.stats:
            prstats("Combination automaton", lcfa)
        if seq is not None:
            with open(seq, "w") as f:
                print(f">>> Writing {seq} in format dot\n")
                lcfa.todot(f)
        if poly is not None:
            reb = rebaser(self.P, vector)
            pofa = lcfa.maptrans(lambda a, v: (a, poly2str(Poly(reb.topoly(v)))))
            with open(poly, "w") as f:
                print(f">>> Writing {poly} in format dot\n")
                pofa.todot(f)
        if cut is None:
            if self.r is None:
                raise ValueError("the substitution is not (X^k-)Pisot!")
            cut = compcut(lcfa, self.P, vector, lineq[1], self.verbose)
        elif self.verbose:
            print(f"Bypassing bound computation! cut={cut}")
        fa = flatten(lcfa, self.P, cut, self.verbose)
        if self.stats:
            prstats(f"Linear equation automaton {lineq}", fa)
        fa = fa.trim()
        if self.stats:
            prstats(f"Trimmed automaton {lineq}", fa)
        if self.verbose:
            obscut = realcut(fa, companion(self.P))
            print(f"Real bound is {obscut} (vs {cut}).")
        if minimize:
            fa = fa.minimize()
            if self.stats:
                prstats("After minimization", fa)
        if dorevector:
            revector(fa, lcfa, self.P)
        if not msd:
            fa = nfa.rev(fa, fa.ma.get(tuple([0] * len(lineq[0]))))
            if self.stats:
                prstats("After reversing", fa)
            fa = fa.minimize()
            if self.stats:
                prstats("After second minimization", fa)
        with open(fname, "w") as f:
            print(f">>> Writing {fname} in format {format}\n")
            if format == "Walnut":
                fa.toWalnut(f)
            elif format == "dot":
                fa.todot(f, vformat=None)
            else:
                raise ValueError(f"Invalid format name '{format}'")

    def trylinear(self, lineq, msd, minimize, vector, cut):
        if vector is None:
            vector = [0] * (self.n - 1) + [1]
        fa = lincomb(self.fa, *lineq)
        if self.stats:
            prstats("Combination automaton", fa)
        if cut is None:
            if self.r is None:
                raise ValueError("the substitution is not (X^k-)Pisot!")
            cut = compcut(fa, self.P, vector, lineq[1], self.verbose)
        elif self.verbose:
            print(f"Bypassing bound computation! cut={cut}")
        fa = flatten(fa, self.P, cut, self.verbose)
        if self.stats:
            prstats(f"Linear equation automaton {lineq}", fa)
        fa = fa.trim()
        if self.stats:
            prstats(f"Trimmed automaton {lineq}", fa)
        obscut = realcut(fa, companion(self.P))
        return obscut

    def gencheck(self, fname, msd, name):
        pre = "msd" if msd else "lsd"
        base = f"{pre}_{name}"
        with open(fname, "w") as f:
            print(f">>> Writing {fname} script\n")
            f.write(
                f"""eval add "?{base} x+y=z":
eval total "?{base} Ax Ay Ez x + y = z";
eval injL "?{base} Ax Ay Az Axx x + y = z & xx + y = z => x = xx";
eval injR "?{base} Ax Ay Az Ayy x + y = z & x + yy = z => y = yy";
eval fun "?{base} Ax Ay Az Azz x + y = z & x + y = zz => z = zz";
eval com "?{base} Ax Ay Az x + y = z <=> y + x = z";
eval distrib "?{base} Ax Ay Az (x + y) + z = x + (y + z)";
eval zeroL "?{base} Ax Ay x+0=y <=> x=y";
eval zeroR "?{base} Ax Ay 0+x=y <=> x=y";
eval succL "?{base} Ax Ay Az Eyy Ezz yy=y+1 & zz=z+1 & x+y=z <=> x+yy=x+zz";
eval succR "?{base} Ax Ay Az Eyy Ezz yy=y+1 & zz=z+1 & x+y=z <=> yy+x=zz+x";
eval basics "?{base} Ea Eaa Eaaa a = 1 & aa = 2 & aaa = 3 & a + aa = aaa & a + 4 = 5 & aa + aaa = 5 & aaa + 5 = 8";
"""
            )
