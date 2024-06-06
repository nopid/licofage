def stra(l):
    return f"({', '.join(map(str,l))})"


def poly2str(P, var="X", mult="", decreasing=True):
    l = list(P)
    res = []
    for i, c in enumerate(l):
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
