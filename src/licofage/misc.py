import numpy as np


def stra(l):
    return f"({', '.join(map(str,l))})"


def poly2str(P, var="X", mult="", decreasing=True):
    l = list(P)
    res = []
    for i, c in enumerate(l):
        if abs(c) < 1e-5:
            c = 0
        try:
            cs = f"{c:+g}"
        except:
            if c < 0:
                cs = f"{c}"
            else:
                cs = f"+{c}"
        if cs == "+0":
            continue
        if cs == "+1":
            if i == 0:
                res.append("+1")
            elif i == 1:
                res.append(f"+{var}")
            else:
                res.append(f"+{var}^{i}")
        elif cs == "-1":
            if i == 0:
                res.append("-1")
            elif i == 1:
                res.append(f"-{var}")
            else:
                res.append(f"-{var}^{i}")
        else:
            if i == 0:
                res.append(f"{cs}")
            elif i == 1:
                res.append(f"{cs}{mult}{var}")
            else:
                res.append(f"{cs}{mult}{var}^{i}")
    if decreasing:
        res = res[::-1]
    res = "".join(res)
    if res == "":
        res = "0"
    elif res[0] == "+":
        res = res[1:]
    return res


import re


def str2poly(s, var="X"):
    s = s.replace(" ", "")
    terms = re.findall(r"([+-]?\d*\*?" + var + r"?\^?\d*)", s)
    terms = [term for term in terms if term]
    coeffs = {}
    for term in terms:
        if "X" not in term:
            coef = int(term)
            power = 0
        else:
            parts = term.split(var)
            if parts[0] in ["", "+", "-"]:
                coef = int(parts[0] + "1")
            elif parts[0] == "-":
                coef = -1
            else:
                coef = int(parts[0].replace("*", ""))

            if "^" in parts[1]:
                power = int(parts[1][1:])
            else:
                power = 1

        if power in coeffs:
            coeffs[power] += coef
        else:
            coeffs[power] = coef
    max_power = max(coeffs.keys())
    coef_list = [0] * (max_power + 1)
    for power, coef in coeffs.items():
        coef_list[power] = coef
    return np.polynomial.Polynomial(coef_list)
