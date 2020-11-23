#!/Users/remy/Applications/PyCharmProjects/fehrest/venv/bin/python

import argparse
from index import InvertedIndex


def main(ref, k):
    # Build the index
    index = InvertedIndex(ref, k)
    index.build()
    index.persist()


if __name__ == "__main__":
    # Parse command
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", required=True, help="The path to your reference genome in FASTA format")
    parser.add_argument("-k", "--kmer_length", default=30, help="kmer/word length")
    args = parser.parse_args()
    main(args.ref, args.kmer_length)
