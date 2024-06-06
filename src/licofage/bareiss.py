"""Bareiss algorithm for integer matrices
https://gist.github.com/vxgmichel/080e9999a1020711f27cd60b5c2d14de"""


def do_piv(a, b, c, d, e):
    x = a * d - b * c
    q, r = divmod(x, e)
    assert r == 0
    return q


def adjugate(a):
    sign = 1
    prev = piv = 1
    n = len(a)
    am = [
        list(map(int, row)) + [0] * i + [1] + [0] * (n - i - 1)
        for i, row in enumerate(a)
    ]
    for c in range(n):
        if am[c][c] == 0:
            for i in range(c + 1, n):
                if am[i][c] != 0:
                    break
            else:
                return 0, 0, None
            am[c], am[i], sign = am[i], am[c], -sign
        prev, piv = piv, am[c][c]
        for i, row in enumerate(am):
            if i == c:
                continue
            am[i] = [do_piv(piv, y, row[c], x, prev) for x, y in zip(row, am[c])]
    det = sign * piv
    return sign, det, [row[n:] for row in am]
