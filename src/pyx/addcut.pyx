#cython: language_level=3, boundscheck=False, wraparound=False
__all__ = ["flatten"]

import cython
from cpython cimport array
import array
import numpy as np
cimport numpy as cnp

cnp.import_array()

ctypedef cnp.int32_t DTYPE_t

cdef dict msh
cdef int curs
cdef list back

cdef ms(int s, int[:] v):
    global msh, curs, back
    cdef tuple tv = (s, tuple(v))
    cdef int rv
    if tv in msh:
        return msh[tv]
    msh[tv] = curs
    rv = curs 
    back.append((s, v))
    curs += 1
    return rv

@cython.binding(True)
def flatten(int q0, int[:] S, int[:] A, int[:] V, int[:] T, int[:] VF, int[:] Sindex,  cnp.ndarray[DTYPE_t, ndim=2] VL, cnp.ndarray[DTYPE_t, ndim=2] M, cnp.ndarray[DTYPE_t, ndim=1] v0, int[:] cut, int vobj):
    global msh, curs, back
    cdef int ssidx
    cdef int ttidx
    cdef int s
    cdef int[:] vs
    cdef int a 
    cdef cnp.ndarray[DTYPE_t, ndim=1] v
    cdef int t
    cdef cnp.ndarray[DTYPE_t, ndim=1] vt
    cdef int i
    cdef array.array RS = array.array('i')
    cdef array.array RA = array.array('i')
    cdef array.array RV = array.array('i')
    cdef array.array RT = array.array('i')
    cdef array.array RF = array.array('i')
    cdef int m = len(S)
    cdef array.array todo = array.array('i')
    cdef set seen = set()
    msh = dict()
    curs = 0
    back = []
    cdef int rq0 = ms(q0, v0)
    todo.append(rq0)
    seen.add(rq0)
    while todo:
        ssidx = todo.pop()
        (s, vs) = back[ssidx]
        if VF[s] == 1 and vs[0] == vobj:
            RF.append(ssidx)
        i = Sindex[s]
        if i > -1:
            while i < m and S[i] == s:
                a = A[i]
                v = VL[V[i]]
                t = T[i]
                vt = np.dot(M, vs) + v
                if all(abs(vt) <= cut):
                    ttidx = ms(t, vt)
                    RS.append(ssidx)
                    RA.append(a)
                    RV.append(V[i])
                    RT.append(ttidx)
                    if not ttidx in seen:
                        todo.append(ttidx)
                        seen.add(ttidx)
                i += 1
    return (curs, rq0, RS, RA, RV, RT, RF)
