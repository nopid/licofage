class mapper:
    def __init__(self):
        self.cur = 0
        self.h = {}
        self.back = []

    def __getitem__(self, s):
        if s in self.h:
            return self.h[s]
        self.h[s] = self.cur
        self.back.append(s)
        self.cur += 1
        return self.h[s]

    def __len__(self):
        return self.cur

    def __iter__(self):
        for x in self.h:
            yield x

    def get(self, s):
        return self.h[s]

    def rev(self, i):
        return self.back[i]

    def copy(self):
        res = mapper()
        res.cur = self.cur
        res.h = self.h.copy()
        res.back = self.back.copy()
        return res


class swapper:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __getitem__(self, v):
        if v == self.x:
            return self.y
        if v == self.y:
            return self.x
        return v
