#cython: language_level=3
# -*- coding: utf-8 -*-

import os
import pickle
import numpy as np
from collections import OrderedDict
from pysam import FastaFile
from utils import canonicalize


class InvertedIndex:
    def __init__(self, ref_path, k, idx_path=None):
        self.k = int(k)
        self.ref_path = ref_path
        self.fa = FastaFile(ref_path)
        self.idx_path = idx_path
        self.idx_name = None
        self.kmer_ranges = None
        self.genome_offsets = None
        self.add_map = self.make_add_map()

    def make_add_map(self):
        """
        Potential method for handling multiple contigs:
        kmer offset = offset in contig + cumsum of previous contig lengths
        """
        length_map = {ref: self.fa.get_reference_length(ref) for ref in self.fa.references}
        added, add_map = 0, {}
        for ref in self.fa.references:
            add_map[ref] = added
            added += length_map[ref]
        return add_map

    def build(self):
        """Build the inverted index over each contig present in the reference fasta"""
        cdef int k = self.k
        cdef int i
        pairs = []
        for record in self.fa.references:
            # added = self.add_map[record]
            print("Indexing {}".format(record))
            reference_length = self.fa.get_reference_length(record)
            num_kmers = reference_length - k + 1
            for i in range(num_kmers):
                kmer = self.fa.fetch(record, i, i+k)
                if 'N' in kmer:
                    continue # This will save space and avoid spurious alignments downstream
                canon_kmer = canonicalize(kmer).encode()
                #pairs.append((kmer, i+added))
                pairs.append((canon_kmer, i))
        print("Sorting kmers...")
        pairs.sort()
        kmers, offsets = zip(*pairs)
        self.genome_offsets = np.array(offsets, dtype='uint32')
        del offsets
        self.kmer_ranges = OrderedDict.fromkeys(kmers)
        print("Indexing kmer Offsets")
        cur = kmers[0]
        for i, kmer in enumerate(kmers):
            if cur != kmer:
                self.kmer_ranges[cur] = (start, i)
                cur = kmer
                start = i
        del kmers
        self.persist()

    def query(self, kmer):
        start, finish = self.get_range(kmer)
        return self.genome_offsets[start:finish]

    def get_range(self, kmer):
        if kmer in self.kmer_ranges:
            return self.kmer_ranges[kmer]

    def prepare_disk(self):
        basename = os.path.basename(self.ref_path)
        disk_name = os.path.splitext(basename)[0]+"_idx"
        os.mkdir(disk_name)
        self.idx_name = disk_name
        meta = [self.ref_path, self.k, self.idx_name]
        pickle.dump(meta, open(disk_name+"/"+"meta_info", 'wb'))

    def persist(self):
        basename = os.path.basename(self.ref_path)
        disk_name = self.idx_name+"/kmer_ranges.pkl"
        pickle.dump(self.kmer_ranges, open(disk_name, 'wb'))
        disk_name = self.idx_name+"/genome_offsets.pkl"
        pickle.dump(self.genome_offsets, open(disk_name, 'wb'))
        print("Index for {} written to disk".format(basename))

    def load_indexes(self):
        self.kmer_ranges = pickle.load(open(self.idx_path + '/kmer_ranges.pkl', "rb"))
        self.genome_offsets = pickle.load(open(self.idx_path + '/genome_offsets.pkl', "rb"))
