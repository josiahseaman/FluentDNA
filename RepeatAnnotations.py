"""This module is designed to read in RepeatMasker annotation CSV downloaded in raw table schema
from UCSC.  This annotation contains alignment to the repeat consensus.

It is similar to Annotations.py.  But at the moment, they're separate files because their target
input and output is significantly different."""

# Read annotation file and just mark where things land on the consensus
import sys

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

    def is_good(self):
        return self.rep_end != self.rep_start


def read_repeatmasker_csv(annotation_filename, white_list_key=None, white_list_value=None, strict=True):
    # example: chr20	1403353	1403432	63040735	-	L3	LINE	CR1	3913	3992
    with open(annotation_filename) as csvfile:
        headers = csvfile.readline()  # ensure the header is what I am expecting and matches RepeatAnnotation definition
        assert headers == '#genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd\n', headers
        if white_list_key:
            headers = headers[1:-1].split('\t')  # remove '#' and \n
            white_list_key = headers.index(white_list_key)  # column index matching the filter criteria
        entries = []
        bad_count = 0
        for row in csvfile:
            columns = row.split('\t')
            if white_list_key is not None:
                if not strict and white_list_value in columns[white_list_key] \
                        or strict and white_list_value == columns[white_list_key]:  # [7] = rep_family, [5] = repName
                    annotation = RepeatAnnotation(*columns)
                    if annotation.is_good():
                        entries.append(annotation)
                    else:
                        bad_count += 1
            else:
                entries.append(RepeatAnnotation(*columns))
        print("Discarded", bad_count, "bad entries.")
        assert len(entries), "No matches found for " + str(white_list_value)
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
            line = blank_line_array(consensus_width, 'A')
            for fragment in text_line:
                nucleotides = fragment.genome_span().sample(seq)
                if fragment.strand == '-':
                    nucleotides = rev_comp(nucleotides)
                nucleotides = nucleotides.replace('A', 'Z')  # TEMP: orange color for Skittle at the moment
                if fragment.rep_end < len(nucleotides):  # sequence I have sampled starts before the beginning of the frame
                    nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
                line = line[:fragment.rep_end - len(nucleotides)] + array('u', nucleotides) + line[fragment.rep_end:]
            assert len(line) == consensus_width + 1, display_lines.index(text_line)  # len(line)
            out.write(''.join(line))


def blank_line_array(consensus_width, filler, newline=True):
    return array('u', (filler * consensus_width) + ('\n' if newline else ''))


def write_consensus_sandpile(anno_entries, out_filename, seq):
    consensus_width = max_consensus_width(anno_entries)
    with open(out_filename + '_%i.fa' % consensus_width, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        depth_graph = [0] * consensus_width
        image = [blank_line_array(consensus_width, 'A') for _ in range(len(anno_entries))]
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
    return consensus_width


def unnormalized_histogram_of_breakpoints(anno_entries, out_filename):
    consensus_width = max_consensus_width(anno_entries)
    with open(out_filename + '_%i.fa' % consensus_width, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        depth_graph = [0 for _ in range(consensus_width)]
        image = [array('u', ('A' * consensus_width) + '\n') for _ in range(len(anno_entries))]
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


def histogram_of_breakpoints(anno_entries, out_filename, reference_points=None):
    consensus_width = max_consensus_width(anno_entries)
    with open(out_filename + '_%i.fa' % consensus_width, 'w') as out:
        out.write('>' + just_the_name(out_filename) + '\n')
        breakpoints = [0 for _ in range(consensus_width)]
        coverage = [0 for _ in range(consensus_width)]
        for fragment in anno_entries:
            for x in range(fragment.rep_start, fragment.rep_end + 1):
                coverage[x] += 1  # used for normalization
            x = fragment.rep_start
            # image[depth(breakpoints, x)][x] = 'G'
            breakpoints[x] += 1

            x = fragment.rep_end
            # image[depth(breakpoints, x)][x] = 'C'
            breakpoints[x] += 1

        coverage_to_break_ratio = sum(coverage) / sum(breakpoints)
        # expected = coverage[x] / sum(coverage) * sum(breakpoints)
        expected_breaks = [coverage[x] / coverage_to_break_ratio for x in range(consensus_width)]

        normalized_break_counts = []
        for x in range(consensus_width):
            normalized_break_counts.append(breakpoints[x] / max(1, expected_breaks[x]))  # 566 copies is a good cutoff for L1
        int_break_counts = [min(1000, int(normalized_break_counts[x] * 20.0)) for x in range(consensus_width)]
        greatest_depth = max(int_break_counts) + 1
        image = [array('u', ('A' * consensus_width) + '\n') for _ in range(greatest_depth)]
        for line_number, line in enumerate(image[:greatest_depth]):
            for x in range(consensus_width):
                if reference_points is not None and x in reference_points:
                    image[line_number][x] = 'T'  # can be overwritten by 'G'
                if int_break_counts[x] >= line_number:
                    image[line_number][x] = 'G'
            out.write(''.join(line))


def lengths_by_repeat_name(anno_entries):
    endings = {}
    for entry in anno_entries:
        current = endings[entry.rep_name] if entry.rep_name in endings else 0
        endings[entry.rep_name] = max(current, entry.rep_end)
    return endings


def output_transposon_fasta(ano_entries, out_filename, seq):
    genome, chromosome, family, mode = just_the_name(out_filename).split('_')
    # sort first by repName then by start position inside the same repName
    sorted_entries = [a for a in sorted(ano_entries, key=lambda x: x.rep_name + '{:010d}'.format(x.rep_start))]
    with open(out_filename + '.fa', 'w') as out:
        for fragment in sorted_entries:
            info = "__".join([fragment.rep_name, family, genome, chromosome, str(fragment.geno_start), str(fragment.geno_end), fragment.strand])
            out.write('>' + info + '\n')
            nucleotides = fragment.genome_span().sample(seq)
            if fragment.strand == '-':
                nucleotides = rev_comp(nucleotides)
            out.write(nucleotides + '\n')  # we're currently not dividing long lines (not FASTA standard)


def output_archetypes_fasta(ano_entries, out_filename, seq):
    genome, chromosome, family, mode = just_the_name(out_filename).split('_')
    sorted_entries = [a for a in sorted(ano_entries, key=lambda x: -len(x))]
    archetypes = {}
    for entry in sorted_entries:  # only the first of each rep_name will make it into archetypes
        if entry.rep_name not in archetypes:
            archetypes[entry.rep_name] = entry

    with open(out_filename + '.fa', 'w') as out:
        for fragment in archetypes.values():
            info = "__".join([fragment.rep_name, family, genome, chromosome, str(fragment.geno_start), str(fragment.geno_end), fragment.strand])
            out.write('>' + info + '\n')
            nucleotides = fragment.genome_span().sample(seq)
            if fragment.strand == '-':
                nucleotides = rev_comp(nucleotides)
            out.write(nucleotides + '\n')  # we're currently not dividing long lines (not FASTA standard)
        print("Wrote", out_filename)


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
        ending_points = lengths_by_repeat_name(anno_entries)
        print(ending_points)
        histogram_of_breakpoints(anno_entries, filename, reference_points=set(ending_points.values()))
    elif key == 'raw_breaks':
        unnormalized_histogram_of_breakpoints(anno_entries, filename)
    elif key == 'fasta':
        output_transposon_fasta(anno_entries, filename, seq)
    elif key == 'archetypes':
        output_archetypes_fasta(anno_entries, filename, seq)
    else:
        print("Key not recognized:", key)



def test_reader():
    # Test Reader
    entries = read_repeatmasker_csv(r'data\RepeatMasker_chr20_alignment.csv', 'repFamily', 'CR1')
    # print(str(entries))
    assert str(entries) == open('data\L3_test.txt', 'r').read(), "String representation doesn't match expected.  Did you read in data\RRepeatMasker_chr20_alignment.csv?"


def filter_repeats_by_chromosome(repeat_entries, contig_name):
    return [x for x in repeat_entries if x.geno_name == contig_name]


def filter_repeats_by_chromosome_and_family(repeat_entries, contig_name, family):
    return [x for x in repeat_entries if x.geno_name == contig_name and x.rep_family == family]


if __name__ == '__main__':
    # test_reader()
    annotation_file = r'data\RepeatMasker_all_alignment.csv'  # RepeatMasker_all_alignment.csv'  RepeatMasker_chr20_alignment
    layout_mode, column, rep_name = sys.argv[1], sys.argv[2], sys.argv[3]  # 'condense', 'repName', 'L1PA3'
    # column, rep_name = 'repName', 'L1PA3'  # ( repName 'repFamily', 'ERV1')  # 'TcMar-Tigger, TcMar-Mariner  # 'ERVK, ERV1, ERVL, L1, Alu, MIR
    # mode = 'condense'  # 'breaks' raw_breaks
    rep_entries = read_repeatmasker_csv(annotation_file, column, rep_name, strict=True)
    target_chr = 'chr1' if len(sys.argv) < 5 else sys.argv[4]
    if 'breaks' not in layout_mode:
        sequence = pluck_contig(target_chr, 'data/hg38.fa')
        rep_entries = filter_repeats_by_chromosome(rep_entries, target_chr)
    else:
        sequence = ''
        target_chr = ''
    print("Found %i entries under %s" % (len(rep_entries), str(rep_name)))

    layout_repeats(rep_entries, 'data/hg38_' + target_chr + '_' + rep_name.replace('_', ''), sequence, layout_mode)
    print('Done')
