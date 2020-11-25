# Fehrest: A Seed-and-Extend aligner in Python

![Scripted Name](img/fehrest.png)

Index-assisted alignment of NGS reads to a reference genome is one of my favorite problems in computational biology. While this problem has already been studied at great length, this aligner is a pet project for me to become more familiar with writing high-performance code in python.

The aligner is an extension of an assignment from one of my favorite classes during my degree and follows the classic seed-and-extend search strategy using an inverted index to solve the offline, inexact-matching problem formally stated as follows:

> For a given text, T, and a query string, P, find *all* occurrences of P in T with up to k edits.

Where an edit may be an insertion, deletion, or a mismatch.



## TODO:

Tier 1:
* Extend to handle size of human genome (write postings list to disk)
* Forward & reverse strand processing

Tier 2:

- Unit Testing
- Make pip installable

Tier 3:
* Implement minimizer index (instead of kmers) with seed, chain, and extend approach
* Semi-global alignment (parasail)
* Process multiple seeds in parallel
* Bit encoding of kmers



## Things I've tried:
### Indexing
* Using SeqIO no default dict

  * ```bash
    python index.py -r input/chr22.fasta  481.03s user 100.30s system 92% cpu 10:30.32 total
    ```

* Using numpy arrays
      So slow I just stopped it. Numerous sources seemed to indicate that the array was best

* Using python arrays and default dict

  * ```bash
    python index.py -r input/chr22.fasta  480.56s user 86.01s system 93% cpu 10:05.86 total
    ```

* Using pysam instead of SeqIO: https://pysam.readthedocs.io/en/latest/api.html#fasta-files

  * ```bash
    python index.py -r input/chr22.fasta  174.80s user 53.40s system 93% cpu 4:03.46 total
    ```
  

- Cython!

  - ```bash
    ./idx.py -r input/chr22.fasta -k 5  67.58s user 0.46s system 99% cpu 1:08.64 total
    ```

## Command Tracking

```bash
remy@Remys-MacBook-Pro-2 input % wgsim -N100 -1150 -S12 -h sars_cov_2.fasta sars_cov_2.fastq tmp.fastq
```
