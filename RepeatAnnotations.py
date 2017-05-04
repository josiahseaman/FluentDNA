"""This module is designed to read in RepeatMasker annotation CSV downloaded in raw table schema
from UCSC.  This annotation contains alignment to the repeat consensus.

It is similar to Annotations.py.  But at the moment, they're separate files because their target
input and output is significantly different."""

# Read annotation file and just mark where things land on the consensus

from DDVUtils import pluck_contig, just_the_name, rev_comp
from Span import Span


class RepeatAnnotation(object):
    def __init__(self, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd):
        self.geno_name = genoName
        self.geno_start = int(genoStart)
        self.geno_end = int(genoEnd)
        self.geno_left = int(genoLeft)
        self.strand = strand
        self.rep_name = repName
        self.rep_class = repClass
        self.rep_family = repFamily
        self.rep_start = int(repStart)
        self.rep_end = int(repEnd)

    def __repr__(self):
        return ' '.join([str(x) for x in (self.geno_name, self.geno_start, self.geno_end, self.geno_left, self.strand,
                                          self.rep_name, self.rep_class, self.rep_family, self.rep_start, self.rep_end)])

    def __len__(self):
        return self.geno_end - self.geno_start

    def genome_span(self):
        return Span(self.geno_start, self.geno_end, self.geno_name, self.strand, zero_ok=False)

    def check_length(self):
        geno_size = self.geno_end - self.geno_start
        rep_size = self.rep_end - self.rep_start
        if geno_size != rep_size:
            print(geno_size, rep_size, geno_size - rep_size)


def read_repeatmasker_csv(annotation_filename, white_list=None):
    with open(annotation_filename) as csvfile:
        headers = csvfile.readline()  # ensure the header is what I am expecting and matches RepeatAnnotation definition
        assert headers == '#genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd\n', headers
        entries = []
        for row in csvfile:
            entry = RepeatAnnotation(*row.split('\t'))
            if white_list:
                if entry.rep_name in white_list:
                    entries.append(entry)
            else:
                entries.append(entry)
        return entries


def write_aligned_repeat_consensus(anno_entries, out_filename, seq):
    consensus_width = max(max([e.rep_end for e in anno_entries]),
                          max([e.rep_start + len(e) for e in anno_entries]))  # two different ways of finding the end
    print("Consensus width!", consensus_width)
    with open(out_filename, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        for fragment in anno_entries:
            nucleotides = fragment.genome_span().sample(seq)
            if fragment.strand == '-':
                nucleotides = rev_comp(nucleotides)
            line = 'A' * (fragment.rep_end - len(nucleotides)) + nucleotides
            line += 'A' * (consensus_width - len(line)) + '\n'
            assert len(line) == consensus_width + 1, len(line)
            out.write(line)
    # chr20	1403353	1403432	63040735	-	L3	LINE	CR1	3913	3992


if __name__ == '__main__':
    # Test Reader
    annotation = r'data\RepeatMasker_chr20_alignment.csv'
    entries = read_repeatmasker_csv(annotation, {'L3'})
    assert str(entries) == open('data\L3_test.txt', 'r').read(), "String representation doesn't match expected.  Did you read in data\RepeatMasker_chr12_alignment.csv?"

    # Test Writer
    seq = pluck_contig('chr20', 'data/hg38_chr20.fa')
    write_aligned_repeat_consensus(entries, 'data/hg38_chr20_L3.fa', seq)
