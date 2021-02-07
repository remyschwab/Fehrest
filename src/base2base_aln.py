"""
Implementations of common per-base alignment methods
"""
import numpy as np


def hamming_distance(a, b, d):
    """Return the number of mismatches between two strings of equal length"""
    assert len(a) == len(b)
    mm = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            mm += 1
        if mm > d:
            break
    return mm


def global_align(x, y, m, g):
    """Perform global alignment with:
        m: mismatch penalty
        g: linear gap penalty
    """
    mat = np.zeros((len(x)+1, len(y)+1), dtype=int)
    mat[0, 1:] = range(1, len(y)+1)
    mat[0, 1:] = np.multiply(mat[0, 1:], g)
    mat[1:, 0] = range(1, len(x)+1)
    mat[1:, 0] = np.multiply(mat[1:, 0], g)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = m if x[i-1] != y[j-1] else 0
            if x[i-1] == "." or y[j-1] == ".":
                delt = 0
            mat[i, j] = min(mat[i-1, j-1]+delt, mat[i-1, j]+g, mat[i, j-1]+g)
    return mat[len(x), len(y)]


# TODO: Local Alignment
def local_align():
    pass


# TODO: Semi-Global Alignment
def semi_align():
    pass
