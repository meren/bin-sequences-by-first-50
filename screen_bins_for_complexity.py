# -*- coding: utf-8 -*-

info = """
    This simple script takes the output of bin_based_on_first_50.py, and uses Shannon entropy
    to quantify the uncertainty asssociated with the representative sequences of bins saved
    in that file, so the input file format is expected to look like this:

    24   TAGGACGTAAGCGCACGCAGGCGTTTGTTAAGTCAGATGTGAAATCCCCG
    15   TAGGAAATCTTCGGCAATGGGCGAAAGCCTGACCGAGCAACGCCGCGTGA
    22   TAGGCGTAAAGGGAGCGTAGGCGGATGATTAAGTGGGATGTGAAATACCC
    81   TAGGCGTAAAGGGCACGCAGGCGGTGACTTAAGTGAGGTGTGAAAGCCCC
    19   TAGGCGTAAAGGGTGCGTAGTGGTTTCTTAAGTCAGAGGTGAAAGGCTAC
    31   TAGGCGTAAAGGGTGCGTAGGTGGTTTCTTAAGTCAGAGGTGAAAGGCTA
    11   TAGGCGTAAAGCGCGACGCGAGGTCGGTTTGTTAAGTCAGAGTGTGAAAC
    21   TAGGCGTAAAGGGCACGCAGGCGGTGACTTAAGTGAGGTGTGAAAAGCCC
    10   TAGGCGTAAAGCGAGCGCAGGCGGTTCCTTAAGTCTGATGTGAAAGCCCC
    10   TAGGCGTAAAGGGAGCGTAGGCGGACTTCTAAGTGAGATGTGAAATACCC
    (...)

    In return, this script creates a file with the sequences found, and their
    Shannon entropy value. The smaller the entropy value, the less complex
    the sequence is.

"""

import sys
from scipy import log2 as log


try:
    bins_output = sys.argv[1]
except IndexError:
    print 'An input file is expected. Info about this script:'
    print info
    print 'Exiting.'
    sys.exit(-1)


bins = [line.strip().split()[1] for line in open(sys.argv[1])]


class EntropyError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Error: %s' % self.e


def entropy(l):
    P = lambda n: (len([x for x in l if x.upper() == n.upper()]) * 1.0 / len(l)) + 0.0000000000000000001
    return -(sum([P(N) * log(P(N)) for N in ['A', 'T', 'C', 'G', '-']]))


entropy_bin_tuples = []

def main(fasta_path):
    for i in range(0, len(bins)):
        sys.stderr.write('\rComputing entropy for bins: %s' \
             % (['|', '/', '-', '\\'][i % 4]))
        sys.stderr.flush() 
        entropy_bin_tuples.append((entropy(bins[i]), bins[i]),)
    sys.stderr.write('\n')

    output_file_path = sys.argv[1] + '.entropy'
    output = open(output_file_path, 'w')
    for tpl in sorted(entropy_bin_tuples, reverse = True):
        output.write('%.6f\t%s\n' % (tpl[0], tpl[1]))
    output.close()

    sys.stderr.write("Output file is ready: %s\n" % output_file_path)

if __name__ == '__main__':
    main(sys.argv[1])
