"""Utility functions"""


def bit_encode(nuc):
    pass


def pack_bkmer(seq):
    pass


def canonicalize(seq):
    return min(seq, reverse_complement(seq))


def reverse_complement(seq):
    tab = str.maketrans("ACTG", "TGAC")
    return seq.translate(tab)[::-1]
