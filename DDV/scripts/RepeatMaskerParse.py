"""This module is designed to read in RepeatMasker annotation CSV downloaded in raw table schema
from UCSC.  This annotation contains alignment to the repeat consensus.

It is similar to Annotations.py.  But at the moment, they're separate files because their target
input and output is significantly different.

The purpose of this script is to show
a variety of MSA pulled from RepeatMasker annotations across the genome.
The intent was to visualize the diversity and abundance at each site.
All the reusable functionality of TransposonLayout has been moved to MultipleAlignmentLayout.
A fairly small script could process the repeatMasker annotation file into a folder of fasta MSA
that could be visualized with MultipleAlignmentLayout.
The crucial values are the within repeat coordinates rep_end.  I found experimentally
that rooting the MSA on the last nucleotide of each line (not the start) gave the most
coherent MSA."""

# Read annotation file and just mark where things land on the consensus
import os
import re
import sys
from typing import List
import unicodedata

from DNASkittleUtils.DDVUtils import rev_comp, editable_str
from DNASkittleUtils.Contigs import pluck_contig, just_the_name, Contig, write_contigs_to_file
from DDV import gap_char
from DDV.Span import Span


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
        if isinstance(white_list_value, list):
            strict = False  # the 'in' operator will check for exact matches inside the list
        if white_list_value is not None:
            headers = headers[1:-1].split('\t')  # remove '#' and \n
            white_list_key = headers.index(white_list_key)  # column index matching the filter criteria
        entries = []
        bad_count = 0
        for row in csvfile:
            columns = row.split('\t')
            if white_list_value is not None:
                if not strict and white_list_value in columns[white_list_key] \
                        or strict and white_list_value == columns[white_list_key]:  # [7] = rep_family, [5] = repName
                    bad_count = add_good_annotation(bad_count, columns, entries)
            else:
                bad_count = add_good_annotation(bad_count, columns, entries)
        print("Discarded", bad_count, "bad entries.")
        assert len(entries), "No matches found for " + str(white_list_value)
        return entries


def add_good_annotation(bad_count, columns, entries):
    annotation = RepeatAnnotation(*columns)
    if annotation.is_good():
        entries.append(annotation)
    else:
        bad_count += 1
    return bad_count


def grab_aligned_repeat(consensus_width, contig, fragment):
    line = blank_line_array(consensus_width, gap_char, newline=False)
    nucleotides = fragment.genome_span().sample(contig.seq)
    if fragment.strand == '-':
        nucleotides = rev_comp(nucleotides)
    if fragment.rep_end - len(nucleotides) < 0:  # sequence I have sampled starts before the beginning of the frame
        nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
    line = line[:fragment.rep_end - len(nucleotides)] + editable_str(nucleotides) + line[fragment.rep_end:]
    assert len(line) == consensus_width, "%i, %i" % (len(line), consensus_width, )

    return line


def blank_line_array(consensus_width, filler, newline=True):
    return editable_str((filler * consensus_width) + ('\n' if newline else ''))


def max_consensus_width(anno_entries):
    consensus_width = max(max([e.rep_end for e in anno_entries]),
                          max([e.rep_start + len(e) for e in anno_entries]))  # two different ways of finding the end
    return consensus_width


def valid_filename(value):
    """ Normalizes string, converts to lowercase, removes non-alpha characters,
    and converts spaces to hyphens. """
    value = unicodedata.normalize('NFKD', value)
    value = re.sub('[^\w\s-]', '', value).strip()
    value = re.sub('[-\s]+', '-', value)
    return value

def output_transposon_folder(anno_entries: List[RepeatAnnotation], chr_contigs, out_folder_name):
    os.makedirs(out_folder_name, exist_ok=True)
    # Make each chr accessible
    chr_dict = {c.name: c for c in chr_contigs}
    #Each repName gets one file
    rep_names = list({x.rep_name for x in anno_entries})
    for rep_name in rep_names:
        annotations = [x for x in anno_entries if x.rep_name == rep_name] #filter
        annotations.sort(key=lambda x: x.rep_start)

        consensus_width = max_consensus_width(annotations)
        all_copies_of_this_type = []
        for rep in annotations:
            processed_seq = grab_aligned_repeat(consensus_width, chr_dict[rep.geno_name], rep)
            contig_name = '%s__%s:%i-%i' %(rep.rep_name, rep.geno_name, rep.geno_start, rep.geno_end)
            c = Contig(contig_name, ''.join(processed_seq))
            all_copies_of_this_type.append(c)
        if len(all_copies_of_this_type) > 1:
            rep_filename = '__'.join(
                [annotations[0].rep_name, annotations[0].rep_family, annotations[0].rep_class])
            write_contigs_to_file(os.path.join(out_folder_name,
                                               valid_filename(rep_filename) + '.fa'),
                                  all_copies_of_this_type,
                                  verbose=True)
        else:
            print("Skipping", rep_name, "because there was only one copy.")


def filter_repeats_by_chromosome(repeat_entries, contig_name):
    return [x for x in repeat_entries if x.geno_name == contig_name]


def filter_repeats_by_chromosome_and_family(repeat_entries, contig_name, family):
    return [x for x in repeat_entries if x.geno_name == contig_name and x.rep_family == family]

def filter_simple_repeats(repeat_entries, return_only_simple_repeats=False):
    before = len(repeat_entries)
    if return_only_simple_repeats:
        repeat_entries = [x for x in repeat_entries if x.rep_class == 'Simple_repeat']  # keep only 'simple'
    else:
        repeat_entries = [x for x in repeat_entries if x.rep_class != 'Simple_repeat']  # remove 'simple'
    difference = before - len(repeat_entries)
    print("Removed", difference, "repeats", "{:.1%}".format(difference / before), "of the data.")
    return repeat_entries


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Required arguments: input_repeat_masker_csv input_fasta output_dir chromosome")
        print('Example: supplemental/RepeatMasker_chr20_alignment.csv "/Genomes/Human/hg38.fa" supplemental/Hg38_chr20_repeats chr20')
    annotation_file = sys.argv[1]  # "supplemental\RepeatMasker_all_alignment.csv"
    input_fasta = sys.argv[2]
    output_dir = sys.argv[3]  # "supplemental\Hg38_chr18_repeats"
    target_chr = sys.argv[4]  # "chr18"

    column = 'repName'
    # column, rep_name = 'repName', 'L1PA3'  # ( repName 'repFamily', 'ERV1')  # 'TcMar-Tigger, TcMar-Mariner  # 'ERVK, ERV1, ERVL, L1, Alu, MIR
    rep_name = None if len(sys.argv) < 6 else sys.argv[5]  # "L1PA3"
    rep_entries = read_repeatmasker_csv(annotation_file, column, rep_name, strict=True)
    sequence = pluck_contig(target_chr, input_fasta)
    rep_entries = filter_repeats_by_chromosome(rep_entries, target_chr)
    print("Found %i entries under %s" % (len(rep_entries), str(rep_name) if rep_name else 'All'))
    rep_entries = filter_simple_repeats(rep_entries, return_only_simple_repeats=False)
    print("Found %i non-simple entries after filtering" % len(rep_entries))

    output_transposon_folder(rep_entries, [Contig(target_chr, sequence)], output_dir)
    print('Done')
