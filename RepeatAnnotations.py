"""This module is designed to read in RepeatMasker annotation CSV downloaded in raw table schema
from UCSC.  This annotation contains alignment to the repeat consensus.

It is similar to Annotations.py.  But at the moment, they're separate files because their target
input and output is significantly different."""

# Read annotation file and just mark where things land on the consensus

import csv


class RepeatAnnotation(object):
    def __init__(self, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd):
        self.genoName = genoName
        self.genoStart = int(genoStart)
        self.genoEnd = int(genoEnd)
        self.genoLeft = int(genoLeft)
        self.strand = strand
        self.repName = repName
        self.repClass = repClass
        self.repFamily = repFamily
        self.repStart = int(repStart)
        self.repEnd = int(repEnd)

    def __repr__(self):
        return ' '.join([str(x) for x in (self.genoName, self.genoStart, self.genoEnd, self.genoLeft, self.strand,
                        self.repName, self.repClass, self.repFamily, self.repStart, self.repEnd)])

    def __len__(self):
        return self.repEnd - self.repStart


def read_repeatmasker_csv(annotation_filename, white_list=None):
    with open(annotation_filename) as csvfile:
        headers = csvfile.readline()  # ensure the header is what I am expecting and matches RepeatAnnotation definition
        assert headers == '#genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd\n', headers
        entries = []
        for row in csvfile:
            entry = RepeatAnnotation(*row.split('\t'))
            if white_list:
                if entry.repName in white_list:
                    entries.append(entry)
            else:
                entries.append(entry)
        return entries


def write_aligned_repeat_consensus(entries, out_filename):
    consensus_width = max([e.repEnd for e in entries])
    with open(out_filename, 'w') as out:
        out.write('>' + out_filename + '\n')
        for fragment in entries:
            line = 'X' * fragment.repStart + 'A' * len(fragment)
            line += 'X' * (consensus_width - len(line)) + '\n'
            assert len(line) == consensus_width + 1, len(line)
            out.write(line)

if __name__ == '__main__':
    # Test Reader
    annotation = r'data\RepeatMasker_chr12_alignment.csv'
    entries = read_repeatmasker_csv(annotation, {'L3'})
    assert str(entries) == open('data\L3_test.txt', 'r').read(), "String representation doesn't match expected.  Did you read in data\RepeatMasker_chr12_alignment.csv?"

    # Test Writer
    write_aligned_repeat_consensus(entries, 'hg38_chr12_L3.fa')