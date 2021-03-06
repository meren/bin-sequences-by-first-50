In this small repository there are two scripts:

* `bin_based_on_first_50.py`
* `screen_bins_for_complexity.py`


## bin_based_on_first_50.py

`bin_based_on_first_50.py` takes a FASTA file path as an argument and creates bins of reads, in respect to the first 50 nucleotides of them, which makes it possible to answer the question _'what is the number of sequences that have the same first 50 nucleotides'_.

If the FASTA file contains shotgun metagenomic reads, the mean number of reads that share the _same_ first 50 nucleotides in each bin would expected to be very small. If you have large bins, it might indicate adapter contamination or other problems.

After obtaining files in this directory, it could be run as follows:

     python bin_based_on_first_50.py /path/to/fasta.fa

By default, the code only analyses the first 10,000,000 reads in the file. Limit can be changed by editing `limit` variable `bin_based_on_first_50.py`.

An example output should look like this:

     meren ~ $ python bin_based_on_first_50.py HMP2102M.fasta
     Approximate number of entries that have been processed so far: ~10,000,000
     Processing bins dictionary: /
     Done.
     
     Results:
     Number of sequences processed                   :  10,000,000
     Number of sequences that were shorter than 50   :  245,279
     Number of sequences that were taken into account:  9,754,721
     Number of bins based on first 50 bases          :  7,817,621
     Number of bins with only one sequence           :  6,656,440 (85.15% of all bins)
     Mean number of sequences in each bin            :  1.24778625109
     Standard deviation                              :  0.835146876261
     Maximum number of sequence in one bin           :  155 (0.001589% of all sequences)
     
     (unique first 50 nucleotides and their counts are stored in ./unique_first_50s)

And this is the output file that was mentioned above:

     meren ~ $ wc -l unique_first_50s 
     7817621 unique_first_50s
     meren ~ $ head -n 10 unique_first_50s 
     1   GCGTTGCGTTCCGTCCATCGGGGGATGGGCTTCTTCAATCTTGAACCTGT
     1   GCGTTGCGTTTCATCATACTCTACATCCTTTGAAACCTGAAACAAAGCAT
     3   GCGTTGCGTTTGCGTCCGTTTCTGGCTTTCCGCAGTGTCCGTCAGTATAC
     3   GCGTTGCGTTAAAGTCCACATACATCTTGCCGTAAACCATAGCCCGCTTG
     1   GCGTTGCGTTCGTCAAACTCGATCACTACCATACCCGATTCCTTAGGCAA
     1   GCGTTAAATGCGTATGGCTTATTGTTAGTCCACGAAGCCACATTCCACGC
     1   GCGTTAAATGGAATCATGTGATCACCTGTATAACCAGTACGCGGTGCTTT
     1   GCGTTAAATGTATTCACTGCCTTAATAAGGCCTATGGTTCCCACGGACAG
     1   GCGTTAAATGCCGTCAAGACAAAACCGAAACTGATGATTTCAGCTTCCGC
     1   GCGTTGCCCGTTGCAATCATGGGCACCACGCAGATGGACTCAATCCTGCT



## screen_bins_for_complexity.py

The second script is `screen_bins_for_complexity.py` which takes the output of the previous script as an input, and quantifies the uncertainty associated with the identified 50 nt bins to give an idea about the complexity of a bin. Due to the nature of the analysis, reads with repetitive regions or with homopolymer characteristics will result in very low entropy. It can be run simply like this:

     meren ~ $ python screen_bins_for_complexity.py ./unique_first_50s
     Computing entropy for bins: \
     Output file is ready: ./unique_first_50s.entropy

And following is an example output file:

     meren ~ $ tail ./unique_first_50s.entropy
     1.987525   TAGGCGTAAAGCGCACGCAGGCGTTTGTTAAGTCAGATGTGAAATCCCCG
     1.987525   TAGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCC
     1.987525   TAGGACGTAAGCGCACGCAGGCGTTTGTTAAGTCAGATGTGAAATCCCCG
     1.981315   TAGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGGTGTGAAATCCCC
     1.846562   TAGGCGTAAAGAGGGAGCAGGCGGCAGCAAGGGTCTGTGGTGAAAGCCTG
     1.409800   AAAAAAAAAAAAAAAAAAAAAAAAAAAGTTAAGTCAGATGTGAAATCCCC
     1.168141   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCAGATGTGAAATCCCC
     0.874310   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGAAATCCCC
     0.541188   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCCCC

You can send your remarks to meren / mbl.edu.
