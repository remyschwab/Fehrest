#cython: language_level=3
# -*- coding: utf-8 -*-

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
        self.dir = None

    def build(self):
        """Build the inverted index over each contig present in the reference fasta"""
        cdef int k = self.k
        cdef int i
        for record in self.ref.references:
            print("Indexing {}".format(record))
            contig_idx = defaultdict(partial(array, 'I'))
            reference_length = self.ref.get_reference_length(record)
            ref_mlen = reference_length / 1000000
            for i in range(reference_length-k+1):
                if i % 1000000 == 0 and i != 0:
                    sys.stdout.write('\rProcessed %0.0f million kmers out of %0.0f' % (float(i)/1000000, ref_mlen))
                    sys.stdout.flush()
                kmer = self.ref.fetch(record, i, i+k).encode()
                if b'N' in kmer:
                    continue # This will save space and avoid spurious alignments downstream
                contig_idx[kmer].append(i)
            self.persist(record, contig_idx)

    def query(self, kmer, contig_idx):
        # contig_idx = self.load_contig(contig)
        if kmer in contig_idx:
            return contig_idx[kmer]

    def prepare_disk(self):
        basename = os.path.basename(self.ref_path)
        disk_name = os.path.splitext(basename)[0]+"_idx"
        os.mkdir(disk_name)
        self.dir = disk_name
        meta = [self.ref_path, self.k]
        pickle.dump(meta, open(disk_name+"/"+"meta_info", 'wb'))

    def persist(self, ref_name, cntg_idx):
        basename = ref_name
        disk_name = self.dir+"/"+basename+".pkl"
        pickle.dump(cntg_idx, open(disk_name, 'wb'))
        print("\nIndex for {} written to disk".format(ref_name))
