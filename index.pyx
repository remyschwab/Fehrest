#cython: language_level=3
# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import sys
import pickle
from array import array
from collections import defaultdict
from functools import partial
from pysam import FastaFile


class InvertedIndex:
    def __init__(self, ref_path, k):
        self.k = int(k)
        self.ref_path = ref_path
        self.ref = FastaFile(ref_path)
        self.index = {}

    def build(self):
        """Build the inverted index over each contig present in the reference fasta"""
        cdef int k = self.k
        cdef int i
        for record in self.ref.references:
            print("Indexing {}".format(record))
            contig_idx = defaultdict(partial(array, 'I'))
            for i in range(self.ref.get_reference_length(record)-k+1):
                if i % 1000000 == 0 and i != 0:
                    print("Processed {} kmers\r".format(i))
                kmer = self.ref.fetch(record, i, i+k).encode()
                if b'N' in kmer:
                    continue # This will avoid spurious alignments downstream
                contig_idx[kmer].append(i)
            print("Processed {} kmers".format(i))
            self.index[record] = contig_idx

    def query(self, kmer, contig):
        contig_idx = self.index[contig]
        if kmer in contig_idx:
            return contig_idx[kmer]

    def persist(self):
        basename = os.path.basename(self.ref_path)
        disk_name = os.path.splitext(basename)[0]+".pkl"
        pickle.dump([self.k, self.index, self.ref_path], open(disk_name, 'wb'))
        print("Index written to disk")

    # TODO
    # def compute_minimzer(self, read, n, k):
    #     Qi = deque()
    #     for i in range(k):
    #         while Qi and arr[i] >= arr[Qi[-1]]:
    #             Qi.pop()
    #         Qi.append(i)
    #     for i in range(k, n):
    #         print(str(arr[Qi[0]]) + " ", end="")
    #         while Qi and Qi[0] <= i - k:
    #             Qi.popleft()
    #         while Qi and arr[i] >= arr[Qi[-1]]:
    #             Qi.pop()
    #         Qi.append(i)
    #     print(str(arr[Qi[0]]))
