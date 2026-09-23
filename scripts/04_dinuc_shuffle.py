#!/usr/bin/env python3
"""
04_dinuc_shuffle.py
-------------------
Altschul-Erickson (1985) dinucleotide-preserving sequence randomisation algorithm.
Generates a random sequence with IDENTICAL mononucleotide and dinucleotide frequencies
to the original input sequence while destroying higher-order sequence patterns.
Maintains the Eulerian-path structure of the nucleotide transition graph.

Reference:
Altschul, S. F., & Erickson, B. W. (1985). Significance of nucleotide sequence alignments:
a method for random sequence permutation that preserves dinucleotide and codon usage.
Mol. Biol. Evol., 2(6), 526–538.
"""

import random
from collections import Counter

def dinuc_shuffle(seq, rng=None):
    """Generates a randomized sequence preserving mononucleotide and dinucleotide frequencies."""
    rng = rng or random
    seq = seq.upper()
    n = len(seq)
    if n < 3:
        return seq

    # Collect outgoing edges
    edges = {}
    for i in range(n - 1):
        a, b = seq[i], seq[i + 1]
        edges.setdefault(a, []).append(b)

    # Preserve last outgoing edge of each node to guarantee valid Eulerian path
    last_edge = {}
    for c, lst in edges.items():
        last_edge[c] = lst[-1]
        head = lst[:-1]
        rng.shuffle(head)
        edges[c] = head + [lst[-1]]

    # Walk Eulerian path
    pos = {c: 0 for c in edges}
    result = [seq[0]]
    cur = seq[0]
    for _ in range(n - 2):
        nxt = edges[cur][pos[cur]]
        pos[cur] += 1
        result.append(nxt)
        cur = nxt
    result.append(last_edge[cur])
    return "".join(result)

def self_test():
    """Verify algorithm invariants: length, composition, dinucleotides, endpoints."""
    def dinuc_counts(s):
        return Counter(s[i:i+2] for i in range(len(s)-1))

    r = random.Random(42)
    test_seqs = [
        "ATGCATGCATGCGGCCATGATGCATGCATGCGGCC",
        "GGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTT",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "ATGC" * 100,
    ]
    all_ok = True
    for s in test_seqs:
        for _ in range(50):
            shuf = dinuc_shuffle(s, rng=r)
            if (len(shuf) != len(s) or 
                Counter(shuf) != Counter(s) or 
                dinuc_counts(shuf) != dinuc_counts(s) or 
                shuf[0] != s[0] or 
                shuf[-1] != s[-1]):
                all_ok = False
                break
    if all_ok:
        print("All dinucleotide shuffle verification tests PASSED.")
    else:
        print("ERROR: Verification tests failed.")

if __name__ == "__main__":
    self_test()
