"""This module is designed to read in RepeatMasker annotation CSV downloaded in raw table schema
from UCSC.  This annotation contains alignment to the repeat consensus.

It is similar to Annotations.py.  But at the moment, they're separate files because their target
input and output is significantly different."""

# Read annotation file and just mark where things land on the consensus

from DDVUtils import pluck_contig, just_the_name, rev_comp
from Span import Span
from array import array
import math


def int_log(num):
    return int(math.log(num + 1, 1.1))


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
        self.__consensus_span = None
        self.__genome_span = None


    def __repr__(self):
        return ' '.join([str(x) for x in (self.geno_name, self.geno_start, self.geno_end, self.geno_left, self.strand,
                                          self.rep_name, self.rep_class, self.rep_family, self.rep_start, self.rep_end)])

    def __len__(self):
        return self.geno_end - self.geno_start

    def consensus_span(self):
        if self.__consensus_span is None:  # cache object
            self.__consensus_span = Span(self.rep_start, self.rep_end, self.rep_name, '+', zero_ok=False)  # always on + strand
        return self.__consensus_span

    def genome_span(self):
        if self.__genome_span is None:  # cache object
            self.__genome_span = Span(self.geno_start, self.geno_end, self.geno_name, self.strand, zero_ok=False)
        return self.__genome_span

    def check_length(self):
        geno_size = self.geno_end - self.geno_start
        rep_size = self.rep_end - self.rep_start
        if geno_size != rep_size:
            print(geno_size, rep_size, geno_size - rep_size)


def read_repeatmasker_csv(annotation_filename, white_list=None):
    # example: chr20	1403353	1403432	63040735	-	L3	LINE	CR1	3913	3992
    with open(annotation_filename) as csvfile:
        headers = csvfile.readline()  # ensure the header is what I am expecting and matches RepeatAnnotation definition
        assert headers == '#genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd\n', headers
        entries = []
        for row in csvfile:
            entry = RepeatAnnotation(*row.split('\t'))
            if white_list:
                if entry.rep_family in white_list:
                    entries.append(entry)
            else:
                entries.append(entry)
        assert len(entries), "No matches found for " + str(white_list)
        return entries


def condense_fragments_to_lines(anno_entries, crowded_count=10):
    """
    :param anno_entries: list of RepeatAnnotation objects
    :param crowded_count:  point at which you should give up trying to cram more into this line
    :return: list of lists of RepeatAnnotation objects that can fit on one line without overlapping
    """
    if crowded_count < 2:
        return [[entry] for entry in anno_entries]

    lines = [[]]
    for entry in anno_entries:  # type = RepeatAnnotation
        for row, candidate_line in enumerate(lines):
            if len(candidate_line) < crowded_count and all([not entry.consensus_span().overlaps(x.consensus_span()) for x in candidate_line]):
                candidate_line.append(entry)
                break
            elif row == len(lines) - 1:
                lines.append([entry])  # add a new line with this entry in it
                break

    return lines


def write_aligned_repeat_consensus(display_lines, out_filename, seq):
    consensus_width = max(max([e.rep_end for line in display_lines for e in line]),
                          max([e.rep_start + len(e) for line in display_lines for e in line]))  # two different ways of finding the end
    with open(out_filename + '_%i.fa' % consensus_width, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        for text_line in display_lines:
            line = array('u', ('A' * consensus_width) + '\n')
            for fragment in text_line:
                nucleotides = fragment.genome_span().sample(seq)
                if fragment.strand == '-':
                    nucleotides = rev_comp(nucleotides)
                nucleotides = nucleotides.replace('A', 'Z')  # TEMP: orange color for Skittle at the moment
                if fragment.rep_end - len(nucleotides) < 0:  # sequence I have sampled starts before the beginning of the frame
                    nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
                line = line[:fragment.rep_end - len(nucleotides)] + array('u', nucleotides) + line[fragment.rep_end:]
            assert len(line) == consensus_width + 1, display_lines.index(text_line)  # len(line)
            out.write(''.join(line))


def write_consensus_sandpile(anno_entries, out_filename, seq):
    consensus_width = max_consensus_width(anno_entries)
    with open(out_filename + '_%i.fa' % consensus_width, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        depth_graph = [0 for i in range(consensus_width)]
        image = [array('u', ('A' * consensus_width) + '\n') for i in range(len(anno_entries))]
        for fragment in anno_entries:
            nucleotides = fragment.genome_span().sample(seq)
            if fragment.strand == '-':
                nucleotides = rev_comp(nucleotides)
            nucleotides = nucleotides.replace('A', 'Z')  # TEMP: orange color for Skittle at the moment
            if fragment.rep_end - len(nucleotides) < 0:  # sequence I have sampled starts before the beginning of the frame
                nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
            for x_start, c in enumerate(nucleotides):
                x = x_start + (fragment.rep_end - len(nucleotides))
                image[depth_graph[x]][x] = c
                depth_graph[x] += 1
        greatest_depth = max(depth_graph)
        for line in image[:greatest_depth]:
            out.write(''.join(line))


def max_consensus_width(anno_entries):
    consensus_width = max(max([e.rep_end for e in anno_entries]),
                          max([e.rep_start + len(e) for e in anno_entries]))  # two different ways of finding the end
    consensus_width = (consensus_width + (12 - consensus_width % 12))  # make it a multiple of 12
    return consensus_width


def histogram_of_breakpoints(anno_entries, out_filename):
    consensus_width = max_consensus_width(anno_entries)
    with open(out_filename + '_%i.fa' % consensus_width, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        depth_graph = [0 for i in range(consensus_width)]
        image = [array('u', ('A' * consensus_width) + '\n') for i in range(len(anno_entries))]
        for fragment in anno_entries:
            x = fragment.rep_start
            image[depth_graph[x]][x] = 'G'
            depth_graph[x] += 1

            x = fragment.rep_end
            image[depth_graph[x]][x] = 'C'
            depth_graph[x] += 1

        greatest_depth = max(depth_graph)
        for line in image[:greatest_depth]:
            out.write(''.join(line))


def layout_repeats(anno_entries, filename, seq, key='condense'):
    filename = filename + '_' + key
    if key == 'condense':
        display_lines = condense_fragments_to_lines(anno_entries, crowded_count=10)
        write_aligned_repeat_consensus(display_lines, filename, seq)
    elif key == 'repStart':
        display_lines = [[a] for a in sorted(anno_entries, key=lambda x: x.rep_start)]
        write_aligned_repeat_consensus(display_lines, filename, seq)
    elif key == 'repEnd':
        display_lines = [[a] for a in sorted(anno_entries, key=lambda x: -x.rep_end)]
        write_aligned_repeat_consensus(display_lines, filename, seq)
    elif key == 'rank':
        display_lines = [[a] for a in sorted(anno_entries, key=lambda x: -len(x))]
        write_aligned_repeat_consensus(display_lines, filename, seq)
    elif key == 'sandpile':
        ranked_entries = [a for a in sorted(anno_entries, key=lambda x: -len(x))]
        write_consensus_sandpile(ranked_entries, filename, seq)
    elif key == 'breaks':
        histogram_of_breakpoints(anno_entries, filename)


def test_reader():
    # Test Reader
    entries = read_repeatmasker_csv(r'data\RepeatMasker_chr20_alignment.csv', ['CR1'])
    # print(str(entries))
    assert str(entries) == open('data\L3_test.txt', 'r').read(), "String representation doesn't match expected.  Did you read in data\RRepeatMasker_chr20_alignment.csv?"


if __name__ == '__main__':
    annotation = r'data\RepeatMasker_all_alignment.csv'
    white_list = ['L1']  # 'TcMar-Tigger, TcMar-Mariner  # 'ERVK, ERV1, ERVL, L1, Alu, MIR
    rep_entries = read_repeatmasker_csv(annotation, white_list)
    print("Found %i entries under %s" % (len(rep_entries), str(white_list)))

    # Test Writer
    mode = 'breaks'
    sequence = pluck_contig('chr20', 'data/hg38_chr20.fa') if mode != 'breaks' else ''
    layout_repeats(rep_entries, 'data/hg38_' + '_'.join(white_list), sequence, mode)
    print('Done')

    test_reader()
