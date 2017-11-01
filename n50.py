#!/usr/bin/env python

#################################
# Copied from "Dina" at https://github.com/dinovski
# on gist https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c
# Dina's code looked like the best way to computer N50 (which is surprisingly not in a library
# anywhere).  I've modified the code to serve as a module inside DDV.

# calculate N50 from fasta file
# N50 = contig length so that half of the contigs are longer and 1/2 of contigs are shorter

import commands
import sys
import os
from itertools import groupby
import numpy

# from Bio import SeqIO

lengths = []

with open(sys.argv[1]) as fasta:
    # parse each sequence by header: groupby(data, key)
    faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

    for record in faiter:
        # join sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        lengths.append(len(seq))

# N50
# the length of the shortest contig such that the sum of contigs of equal
# length or longer is at least 50% of the total length of all contigs

# sort contigs longest>shortest
all_len = sorted(lengths, reverse=True)
csum = numpy.cumsum(all_len)

print
"N: %d" % int(sum(lengths))
n2 = int(sum(lengths) / 2)

# get index for cumsum >= N/2
csumn2 = min(csum[csum >= n2])
ind = numpy.where(csum == csumn2)

n50 = all_len[ind[0]]
print
"N50: %s" % n50

# N90
nx90 = int(sum(lengths) * 0.90)

# index for csumsum >= 0.9*N
csumn90 = min(csum[csum >= nx90])
ind90 = numpy.where(csum == csumn90)

n90 = all_len[ind90[0]]
print
"N90: %s" % n90

# write lengths to file
with open('all_lengths.txt', 'w') as handle:
    handle.write('\n'.join(str(i) for i in lengths))