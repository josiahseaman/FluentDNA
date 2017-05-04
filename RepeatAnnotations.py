"""This module is designed to read in RepeatMasker annotation CSV downloaded in raw table schema
from UCSC.  This annotation contains alignment to the repeat consensus.

It is similar to Annotations.py.  But at the moment, they're separate files because their target
input and output is significantly different."""

# Read annotation file and just mark where things land on the consensus

from DDVUtils import pluck_contig, just_the_name, rev_comp
from Span import Span
from array import array


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
        if self.rep_end == self.rep_start:
            if self.strand == '+':
                self.rep_end += 1
            else:
                self.rep_start += 1

    def __repr__(self):
        return ' '.join([str(x) for x in (self.geno_name, self.geno_start, self.geno_end, self.geno_left, self.strand,
                                          self.rep_name, self.rep_class, self.rep_family, self.rep_start, self.rep_end)])

    def __len__(self):
        return self.geno_end - self.geno_start

    def consensus_span(self):
        return Span(self.rep_start, self.rep_end, self.rep_name, '+', zero_ok=False)  # always on + strand

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


def condense_fragments_to_lines(anno_entries, consensus_width):
    lines = [[]]
    for entry in anno_entries:  # type = RepeatAnnotation
        for row, candidate_line in enumerate(lines):
            if all([not entry.consensus_span().overlaps(x.consensus_span()) for x in candidate_line]):
                candidate_line.append(entry)
                break
            elif row == len(lines) - 1:
                lines.append([entry])  # add a new line with this entry in it
                break

    return lines


def write_aligned_repeat_consensus(anno_entries, out_filename, seq):
    consensus_width = max(max([e.rep_end for e in anno_entries]),
                          max([e.rep_start + len(e) for e in anno_entries]))  # two different ways of finding the end
    print("Consensus width!", consensus_width)
    with open(out_filename, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        display_lines = condense_fragments_to_lines(anno_entries, consensus_width)
        for text_line in display_lines:
            line = array('u', ('A' * consensus_width) + '\n')
            for fragment in text_line:
                if str(fragment) == 'chr20 197571 197817 64246350 + L1PA12_3end LINE L1 1 237':
                    print(fragment)
                nucleotides = fragment.genome_span().sample(seq)
                if fragment.strand == '-':
                    nucleotides = rev_comp(nucleotides)
                nucleotides = nucleotides.replace('A', 'Z')  # TEMP: orange color for Skittle at the moment
                if fragment.rep_end - len(nucleotides) < 0:  # sequence I have sampled starts before the beginning of the frame
                    nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
                line = line[:fragment.rep_end - len(nucleotides)] + array('u', nucleotides) + line[fragment.rep_end:]
            assert len(line) == consensus_width + 1, display_lines.index(text_line)  # len(line)
            out.write(''.join(line))
    # chr20	1403353	1403432	63040735	-	L3	LINE	CR1	3913	3992


def test_reader():
    # Test Reader
    entries = read_repeatmasker_csv(r'data\RepeatMasker_chr20_alignment.csv', ['L3'])
    assert str(entries) == open('data\L3_test.txt', 'r').read(), "String representation doesn't match expected.  Did you read in data\RepeatMasker_chr12_alignment.csv?"


if __name__ == '__main__':
    annotation = r'data\RepeatMasker_chr20_alignment.csv'
    white_list = ['L1PA12_3end']
    entries = read_repeatmasker_csv(annotation, white_list)

    # Test Writer
    seq = pluck_contig('chr20', 'data/hg38_chr20.fa')
    write_aligned_repeat_consensus(entries, 'data/hg38_chr20_' + '_'.join(white_list) +'.fa', seq)
    print('Done')

    test_reader()
