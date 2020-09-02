# sequence-entropy
This is the perl script we use to calculate sequence entropy.
It is put here, so we can reference it in a paper.

It is documented. Try <code>perldoc entropy.pl</code>

It eats an alignment in fasta format and then prints out a long table in a form suitable for feeding to R or gnuplot for plotting.

It has now been obsoleted. Have a look at https://github.com/andrew-torda/seq_compat which is written in go. It is about 50 times faster for calculating sequence entropy, but can also compare two multiple sequence alignments.
