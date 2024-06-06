#cython: language_level=3, boundscheck=False, wraparound=False
"""DFA minimization

   Cython version of the program DFA_minimizer.cc 
   Retrieved from http://users.jyu.fi/~ava/ on 22 November 2023
   Quick and dirty conversion by N. Ollinger, please read
   Copyright and conditions of use below.

   The original program header reads as follows:

/*
  Copyright Antti Valmari 2012, this commented version 2013.
  This program is from

    Antti Valmari: Fast brief practical DFA minimization,
    Information Processing Letters 112 (2012) 213-217

  You may use and adapt the program for scientific purposes at your own risk,
  but you must give credit to the original author and source. Please negotiate
  with me about uses for other purposes.

  If you do not have access to the above-mentioned publication, please see
  A. Valmari, P. Lehtinen: Efficient minimization of DFAs with partial
  transition functions, Symposium on Theoretical Aspects of Computer Science,
  2008, pp. 645-656, http://drops.dagstuhl.de/volltexte/2008/1328/
  That publication explains part of the background. However, this program is
  much further optimized.

  This program inputs a deterministic finite automaton whose transition
  function is not necessarily full, and outputs the minimal automaton
  accepting the same language. The program also contains the removal of
  irrelevant parts of the DFA.

  This program runs in O(n + m log m) time, where n is the number of states
  and m is the number of defined transitions. If the transitions are given in
  the input such that all transitions with the same label are given together,
  then the transitions with another label, and so on, then the lines
  "#include <algorithm>" and "std::sort( C.E, C.E+mm, cmp );" can be removed,
  improving the running time to O(n + m log n). These should be compared to
  the running time of Hopcroft's algorithm, which is O(nk log n), where k is
  the size of the alphabet.

  This program is also fast in practice, I believe.
*/
"""
__all__ = ["minimize"]

from libc.stdlib cimport qsort
from cpython cimport array
import array
import cython

cdef array.array iat = array.array('i', [])
cdef int w
cdef int nn
cdef int mm
cdef int[:] gM
cdef int[:] gW
cdef int[:] gA
cdef int[:] gF
cdef B
cdef C

cdef class partition:
    cdef public int z
    cdef public int[:] E
    cdef public int[:] L
    cdef public int[:] S
    cdef public int[:] F
    cdef public int[:] P

    def __init__(self, int n):
        z = 1 if n>0 else 0
        cdef int[:] E = array.clone(iat, n, zero=True)
        cdef int[:] L = array.clone(iat, n, zero=True)
        cdef int[:] S = array.clone(iat, n, zero=True)
        cdef int[:] F = array.clone(iat, n, zero=True)
        cdef int[:] P = array.clone(iat, n, zero=True)
        for i in range(n):
            E[i] = i
            L[i] = i
            S[i] = 0
        if n > 0:
            F[0] = 0
            P[0] = n
        self.z = z
        self.E = E
        self.L = L
        self.S = S
        self.F = F
        self.P = P

    def mark(self, int e):
        global w, gM, gW
        cdef int[:] E = self.E
        cdef int[:] L = self.L
        cdef int[:] S = self.S
        cdef int[:] F = self.F
        cdef int[:] P = self.P
        cdef int[:] M = gM
        cdef int[:] W = gW
        cdef int s = S[e]
        cdef int i = L[e]
        cdef int j = F[s]+M[s]
        E[i] = E[j]
        L[E[i]] = i
        E[j] = e
        L[e] = j
        if M[s]==0:
            W[w]=s
            w+=1
        M[s]+=1

    def split(self):
        global w, gM, gW
        cdef int[:] E = self.E
        cdef int[:] L = self.L
        cdef int[:] S = self.S
        cdef int[:] F = self.F
        cdef int[:] P = self.P
        cdef int[:] M = gM
        cdef int[:] W = gW
        cdef int z = self.z
        cdef int s = 0
        cdef int j = 0
        cdef int i
        while w>0:
            w-= 1
            s = W[w]
            j = F[s] + M[s]
            if j == P[s]:
                M[s] = 0
                continue
            if M[s] <= P[s]-j:
                F[z] = F[s]
                P[z] = j
                F[s] = j
            else:
                P[z] = P[s]
                F[z] = j
                P[s] = j
            for i in range(F[z],P[z]):
                S[E[i]] = z
            M[s] = 0
            M[z] = 0
            z += 1
        self.z = z


cdef make_adjacent(int[:] K):
    global nn, mm, gF, gA
    cdef int[:] F = gF
    cdef int[:] A = gA
    cdef int q
    cdef int t
    for q in range(nn+1):
        F[q]=0
    for t in range(mm):
        F[K[t]]+=1
    for q in range(nn):
        F[q+1] += F[q]
    for t in range(mm-1,-1,-1):
        F[K[t]] -= 1
        A[F[K[t]]] = t


cdef int rr

cdef reach(int q):
    global B, rr
    cdef int[:] BE = B.E
    cdef int[:] BL = B.L
    cdef int i = BL[q]
    if i >= rr:
        BE[i] = BE[rr]
        BL[BE[i]] = i
        BE[rr] = q
        BL[q] = rr
        rr += 1

cdef rem_unreachable(int[:] T, int[:] L, int[:] H):
    global B, rr, gF, gA, mm
    cdef int[:] BE = B.E
    cdef int[:] BL = B.L
    cdef int[:] BP = B.P
    cdef int[:] F = gF
    cdef int[:] A = gA
    cdef int j = 0
    make_adjacent(T)
    cdef int i=0
    while i<rr:
        for j in range(F[BE[i]], F[BE[i]+1]):
            reach(H[A[j]])
        i+=1
    j = 0
    for t in range(mm):
        if BL[T[t]] < rr:
            H[j] = H[t]
            L[j] = L[t]
            T[j] = T[t]
            j += 1
    mm = j
    BP[0] = rr
    rr = 0

cdef int *cmpL
cdef int cmp_fun(const void *i, const void *j) noexcept nogil:
    cdef int x = cmpL[(<int*>i)[0]]
    cdef int y = cmpL[(<int*>j)[0]]
    if x < y:
        return -1
    if x == y:
        return 0
    return 1

@cython.binding(True)
def minimize(int n, int m, int q0, int[:] T,  int[:] L, int[:] H,  int[:] FF, sorted=False):
    """Minimize the deterministic finite automaton given as input.

The automaton has [n] states and [m] transitions. States and transitions 
labels are encoded as integers. The initial state is [q0]. The accepting
states are listed in the array [FF]. The transitions are listed in the
arrays [T], [L], [H] of length [m].

A transition [(s,a,t) = (T[i],L[i],H[i])] reads as
"from state [s] reading symbol [a] the automaton enters state [t]".

If the transitions are given in an order such that every transitions
with the same label are given together, setting [sorted] to [True] 
avoids a useless sort and improves the overall complexity 
from O(n + m log m) to O(n + m log n).

All arrays are expected to be provided as objects that exposes writable
buffers of integers. It includes:
  - CPython [array.array] of type code ['i'] as [array('i', [1, 2, 3])]
  - NumPy arrays of dtype [np.intc] as [np.array([1, 2, 3], dtype=np.intc)]

The return value is a tuple [(n, m, q0, T, L, H ,FF)] describing the
output minimal automaton. Every array is encoded as an [array.array] of
type code ['i'].

```
from dfamin import minimize
from array import array

print(
    minimize(
        2,
        4,
        0,
        array("i", [0, 0, 1, 1]),
        array("i", [1, 0, 1, 0]),
        array("i", [0, 1, 0, 1]),
        array("i", [0, 1]),
    )
)
```
    """
    global gM, gW, gF, gA, B, C, nn, mm, rr, w, cmpL
    gL = L
    w = 0
    rr = 0
    nn = n
    mm = m
    cdef int[:] A = array.clone(iat, mm, zero=True)
    cdef int[:] F = array.clone(iat, nn + 1, zero=True)
    gA = A
    gF = F
    B = partition(nn)
    cdef int[:] BE = B.E
    cdef int[:] BL = B.L
    cdef int[:] BP = B.P
    cdef int[:] BS = B.S
    cdef int[:] BF = B.F
    cdef int q
    cdef int i
    # Remove states that cannot be reached from the initial state, and from which final states cannot be reached
    reach(q0)
    rem_unreachable(T, L, H)
    for q in FF:
        if BL[q] < BP[0]:
            reach(q)
    ff = rr
    rem_unreachable(H, L, T)
    # Make initial partition
    gM = array.clone(iat, mm + 1, zero=True)
    gW = array.clone(iat, mm + 1, zero=True)
    cdef int[:] M = gM
    cdef int[:] W = gW
    M[0] = ff
    if ff > 0:
        W[w] = 0
        w += 1
        B.split()
    # Make transition partition
    C = partition(mm)
    cdef int[:] CE = C.E
    cdef int[:] CP = C.P
    cdef int[:] CF = C.F
    cdef int[:] CS = C.S
    cdef int[:] CL = C.L
    cdef int z
    cdef int a
    cdef int t
    if mm > 0:
        if not sorted:
            cmpL = &L[0]
            qsort(&CE[0], mm, sizeof(int), &cmp_fun)
        z = 0
        M[0] = 0
        a = L[CE[0]]
        for i in range(mm):
            t = CE[i]
            if L[t] != a:
                a = L[t]
                CP[z] = i
                z += 1
                CF[z] = i
                M[z] = 0
            CS[t] = z
            CL[t] = i
        CP[z] = mm
        C.z = z + 1
    # Split blocks and cords
    make_adjacent(H)
    cdef int b = 1
    cdef int c = 0
    cdef int j
    while c < C.z:
        for i in range(CF[c], CP[c]):
            B.mark(T[CE[i]])
        B.split()
        c += 1
        while b < B.z:
            for i in range(BF[b], BP[b]):
                for j in range(F[BE[i]], F[BE[i] + 1]):
                    C.mark(A[j])
            C.split()
            b += 1
    # Count the numbers of transitions and final states in the result
    cdef int mo = 0
    cdef int fo = 0
    for t in range(mm):
        if BL[T[t]] == BF[BS[T[t]]]:
            mo += 1
    z = B.z
    for b in range(z):
        if BF[b] < ff:
            fo += 1
    resn = B.z
    resm = mo
    resq0 = BS[q0]
    resff = fo
    cdef array.array resT = array.clone(iat, resm, zero=True)
    cdef array.array resL = array.clone(iat, resm, zero=True)
    cdef array.array resH = array.clone(iat, resm, zero=True)
    cdef array.array resFF = array.clone(iat, resff, zero=True)
    cdef int cur = 0
    for t in range(mm):
        if BL[T[t]] == BF[BS[T[t]]]:
            resT[cur]=BS[T[t]]
            resL[cur]=L[t]
            resH[cur]=BS[H[t]]
            cur += 1
    cur = 0
    z = B.z
    for b in range(z):
        if BF[b] < ff:
            resFF[cur] = b
            cur += 1
    return (resn, resm, resq0, resT, resL, resH, resFF)
