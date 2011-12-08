# -*- coding: utf-8 -*-

import sys
import numpy
import random
import fastalib as u

# change this to a larger number to analyze more sequences,
# put 0 to analyze all sequences in a given file
limit = 9999999

bins_dict = {}
binned_sequence_counts = []

input = u.SequenceSource(sys.argv[1])
output = open('unique_first_50s', 'w')

number_of_sequences_that_were_shorter_than_50 = 0

def pp(n):
    """Pretty print function for very big numbers.."""
    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

while input.next():
    if input.pos % 10000 == 0 or input.pos == 1:
        sys.stderr.write('\rApproximate number of entries that have been processed so far: ~%s' % (pp(input.pos)))
        sys.stderr.flush() 
    
    if len(input.seq) < 50:
        number_of_sequences_that_were_shorter_than_50 += 1
        continue

    if limit and input.pos > limit:
        break

    if bins_dict.has_key(input.seq[0:5]):
        if bins_dict[input.seq[0:5]].has_key(input.seq[5:10]):
            if input.seq[0:50] in bins_dict[input.seq[0:5]][input.seq[5:10]]:
                bins_dict[input.seq[0:5]][input.seq[5:10]][input.seq[0:50]] += 1
            else:
                bins_dict[input.seq[0:5]][input.seq[5:10]][input.seq[0:50]] = 1
        else:
            bins_dict[input.seq[0:5]][input.seq[5:10]] = {input.seq[0:50]: 1}
    else:
        bins_dict[input.seq[0:5]] = {input.seq[5:10]: {input.seq[0:50]: 1}}

sys.stderr.write('\n')

for k1 in bins_dict:
    sys.stderr.write('\rProcessing bins dictionary: %s' \
             % (['|', '/', '-', '\\'][int(random.random() * 10) % 4]))
    sys.stderr.flush() 
    for k2 in bins_dict[k1]:
        for unique_fifty in bins_dict[k1][k2]:
            output.write('%d\t%s\n' % (bins_dict[k1][k2][unique_fifty], unique_fifty))
            binned_sequence_counts.append(bins_dict[k1][k2][unique_fifty])

binned_sequence_counts.sort(reverse = True)
total_number_of_sequences_taken_into_account = input.pos - number_of_sequences_that_were_shorter_than_50
maximum_number_of_sequence_in_one_bin = max(binned_sequence_counts)
number_of_bins_with_only_one_sequence = len([t for t in binned_sequence_counts if t == 1])

sys.stderr.write('\nDone.\n\nResults:\n')


print 'Number of sequences processed                   : ', pp(input.pos)
print 'Number of sequences that were shorter than 50   : ', pp(number_of_sequences_that_were_shorter_than_50) 
print 'Number of sequences that were taken into account: ', pp(total_number_of_sequences_taken_into_account)
print 'Number of bins based on first 50 bases          : ', pp(len(binned_sequence_counts))
print 'Number of bins with only one sequence           :  %s (%.2f%% of all bins)'\
                    % (pp(number_of_bins_with_only_one_sequence),
                       number_of_bins_with_only_one_sequence * 100.0 / len(binned_sequence_counts))
print 'Mean number of sequences in each bin            : ', numpy.mean(binned_sequence_counts)
print 'Standard deviation                              : ', numpy.std(binned_sequence_counts)
print 'Maximum number of sequence in one bin           :  %s (%f%% of all sequences)'\
                    % (pp(maximum_number_of_sequence_in_one_bin),
                       maximum_number_of_sequence_in_one_bin * 100.0 / total_number_of_sequences_taken_into_account)
print
print '(unique first 50 nucleotides and their counts are stored in ./unique_first_50s)'
