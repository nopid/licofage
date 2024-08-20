import numpy as np
from numpy.polynomial import Polynomial as Poly
from fractions import Fraction


def PolyFrac(*a):
    "construct fraction polynomial from integer coefficients list [a]"
    return Poly(list(map(lambda x: Fraction(x), a))).trim()


def unfrac(P):
    "approximate fraction coefficients of [P] by floats"
    return Poly(list(map(float, P.coef)))


one = PolyFrac(1)
zero = PolyFrac(0)
x = PolyFrac(0, 1)


def monic(P):
    "Monic poly from [P]"
    return (Fraction(1, 1) / P.coef[-1]) * P


def bm(a):
    """Use the Berlekamp-Massey algorithm to find the minimal polynomial
    of a linearly recurrence sequence of integers a.
    As described in https://arxiv.org/pdf/2211.11721.pdf p. E4
    """
    assert len(a) % 2 == 0, "sequence length has to be even"
    n = len(a) // 2
    R1, R0 = PolyFrac(*a[::-1]), x ** (2 * n)
    if R1 == zero:
        return one
    V1, V0 = one, zero
    while R1.degree() >= n:
        R0, Q, R1 = R1, *divmod(R0, R1)
        V0, V1 = V1, V0 - Q * V1
    return monic(V1)


def naivegcd(P, Q):
    "Compute monic GCD of P and Q"
    a, b = P, Q
    while b.degree() > 0 or b.coef[0] != 0:
        a, b = b, a % monic(b)
    return monic(a)


def naivelcm(P, Q):
    "Compute LCM of P and Q"
    return P * (Q // naivegcd(P, Q))


def commonpoly(*aa):
    "Compute LCM poly of [aa] sequences"
    cur = bm(aa[0])
    for a in aa[1:]:
        nxt = bm(a)
        cur = naivelcm(cur, nxt)
    return cur


def pisot_root(P, prec=1e-4):
    "Compute the Pisot root of [P] if any"
    r = list(unfrac(P).roots())
    r.sort(key=lambda v: abs(v), reverse=True)
    if r[0].real > 1 - prec and abs(r[0].imag) < prec:
        for o in r[1:]:
            if abs(o) >= 1 - prec or abs(o) == 0:
                return None
        return r[0].real
    return None


def gen_pisot_root(P, prec=1e-4):
    "Compute the Pisot root of [P/X^k] if any"
    r = list(unfrac(P).roots())
    r.sort(key=lambda v: abs(v), reverse=True)
    k = 0
    while abs(r[-1]) < prec:
        k += 1
        r.pop(-1)
    if r[0].real > 1 - prec and abs(r[0].imag) < prec:
        for o in r[1:]:
            if abs(o) >= 1 - prec or abs(o) == 0:
                return None
        return k, r[0].real
    return None


def companion(P):
    "Compute companion matrix of integer coefficents [P]"
    d = []
    n = P.degree()
    for x in -P.coef[:-1]:
        assert x.denominator == 1, "not an integer polynomial"
        d.append(x.numerator)
    m = [[0] * i + [1] + [0] * (n - i - 1) for i in range(1, n)]
    return np.array(m + [d])
