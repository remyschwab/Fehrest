#! /Users/remy/Applications/PyCharmProjects/fehrest/venv/bin/python


import os
import argparse
import pickle
from pysam import FastaFile
from array import array
from collections import defaultdict
from functools import partial


class InvertedIndex:
    def __init__(self, ref_path, k):
        self.k = int(k)
        self.ref_path = ref_path
        self.ref = FastaFile(ref_path)
        self.index = {}

    def build(self):
        """Build the inverted index over each contig present in the reference fasta"""
        k = self.k
        for record in self.ref.references:
            print("Indexing {}".format(record))
            contig_idx = defaultdict(partial(array, 'I'))
            for i in range(self.ref.get_reference_length(record)-k+1):
                if i % 1000000 == 0 and i != 0:
                    print("Processed {} kmers".format(i), end="\r")
                kmer = self.ref.fetch(record, i, i+k)
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
    def compute_minimzer(self):
        pass


if __name__ == "__main__":
    # Parse command
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", required=True, help="The path to your reference genome in FASTA format")
    parser.add_argument("-k", "--kmer_length", default=30, help="kmer/word length")
    args = parser.parse_args()
    # Build the index
    index = InvertedIndex(args.ref, args.kmer_length)
    index.build()
    index.persist()
