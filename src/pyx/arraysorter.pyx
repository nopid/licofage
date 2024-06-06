#cython: language_level=3, boundscheck=False, wraparound=False
"""Misc array utility functions."""

__all__ = ["arraysorter", "indexer"]
from libc.stdlib cimport qsort, malloc, free
from libc.string cimport memcpy
from cpython cimport array
import array
import cython

cdef array.array iat = array.array('i', [])

cdef int **cmpdata
cdef int cmpn
cdef int cmp_fun(const void *i, const void *j) noexcept nogil:
    cdef int x = (<int*>i)[0]
    cdef int y = (<int*>j)[0]
    cdef int k = 0
    while k < cmpn and cmpdata[k][x] == cmpdata[k][y]:
        k += 1
    if k == cmpn:
        return 0
    if cmpdata[k][x] < cmpdata[k][y]:
        return -1
    return 1

@cython.binding(True)
def arraysorter(list arraylist):
    "Sort a list of int arrays encoding tuples according to the array values in order."
    global cmpdata, cmpn
    cdef int k = len(arraylist)
    cdef int n = len(arraylist[0])
    cmpn = k
    cdef int[:] ref
    cdef int[:] idx = array.clone(iat, n, zero=False)
    cdef int[:] tmp = array.clone(iat, n, zero=False)
    cdef int *cur
    cdef int i
    cdef int j
    cmpdata = <int **>malloc(k * sizeof(int *))
    for i in range(k):
        ref = arraylist[i]
        cmpdata[i] = &ref[0]
    for i in range(n):
        idx[i] = i
    qsort(&idx[0], n, sizeof(int), &cmp_fun)
    for i in range(k):
        cur = cmpdata[i]
        for j in range(n):
            tmp[j] = cur[idx[j]]
        memcpy(&cur[0], &tmp[0], n*sizeof(int))
    free(cmpdata)


@cython.binding(True)
def indexer(array.array a, int n):
    "Compute first index in array [a] for values [0] to [n]."
    cdef array.array resa = array.clone(iat, n, zero=False)
    cdef int[:]  res = resa
    cdef int cur = -1
    cdef int last = 0
    cdef int m = len(a)
    cdef int[:] data = a
    cdef int i
    cdef int j
    for i in range(m):
        last = cur
        cur = data[i]
        if cur == last:
            continue
        for j in range(last+1, cur):
            res[j] = -1
        res[cur] = i
    return resa
