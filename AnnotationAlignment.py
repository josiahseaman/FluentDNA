import os

from ChainParser import ChainParser
from Span import alignment_chopping_index, AlignedSpans, Span
from DDVUtils import Contig
from RepeatAnnotations import max_consensus_width, read_repeatmasker_csv
from TransposonLayout import TransposonLayout, grab_aligned_repeat, filter_repeats_by_chromosome


# class AnnotationAlignment:
#     def __init__(self):


def create_aligned_annotation_fragments(alignment, repeat_entries):
    aligned_repeats = []  # list of AlignedSpans
    consensus_width = max_consensus_width(repeat_entries)
    discarded = 0
    for repeat in repeat_entries:
        fragment = repeat.genome_span()
        scrutiny_index = max(alignment_chopping_index(alignment, fragment) - 1, 0)  # Binary search
        current = alignment[scrutiny_index]
        # take only the piece that is annotated and add it to aligned_repeats
        if fragment.begin in current.ref:
            if fragment.end in current.ref:  # completely inside aligned area
                # ref[genStart: len]
                ref = Span(fragment.begin, fragment.end, current.ref.contig_name, current.ref.strand)
                # query[queryStart + (genStart - refStart) : len]
                ref_beginning_discarded = fragment.begin - current.ref.begin  # how much of ref was chopped off at the beginning?
                query_start = current.query.begin + ref_beginning_discarded
                query = Span(query_start, query_start + len(ref), current.query.contig_name, current.query.strand)

                # add gaps based on repStart and repEnd
                # TODO: try using rep_end instead  ,repeat.rep_end - len(ref))
                pre_spacer = Span(0, repeat.rep_start, query.contig_name, query.strand)
                post_spacer = Span(0, consensus_width - len(ref) - len(pre_spacer), query.contig_name, query.strand)
                aligned_repeats.append(AlignedSpans(pre_spacer, pre_spacer, 0, 0, is_hidden=True))  # blank space at beginning of line
                aligned_repeats.append(AlignedSpans(ref, query, 0, 0, current.is_master_chain))
                aligned_repeats.append(AlignedSpans(post_spacer, post_spacer, 0, 0, is_hidden=True))  # blank space at end
            else:
                discarded += 1
    print("Discarded %s: %i out of %i partial overlaps" % ("{:.0%}".format(discarded / len(repeat_entries)), discarded, len(repeat_entries)))
    return aligned_repeats


def align_annotation(annotation_filename, ref_fasta, query_fasta, chain_file):
    chromosome = 'chr20'
    # We want to do pieces of the contents of _parse_chromosome_in_chain
    chain = ChainParser(chain_file, ref_fasta, query_fasta, '')
    names, ref_chr = chain.setup_for_reference_chromosome(chromosome)
    chain.create_alignment_from_relevant_chains(ref_chr)

    # modify chain.alignment to only contain annotated stretches
    repeat_entries = read_repeatmasker_csv(annotation_filename, 'repName', 'L1PA3')
    repeat_entries = filter_repeats_by_chromosome(repeat_entries, chromosome)
    consensus_width = max_consensus_width(repeat_entries)
    print("Consensus Width: ", consensus_width)
    chain.alignment = create_aligned_annotation_fragments(chain.alignment, repeat_entries)
    chain.create_fasta_from_composite_alignment()  # populates chain.query_seq_gapped and ref_seq_gapped
    # or if no fasta output is wanted: query_uniq_array, ref_uniq_array = chain.compute_unique_sequence()

    names['ref_gapped'], names['query_gapped'] = chain.write_gapped_fasta(names['ref'], names['query'], prepend_output_folder=True)
    print(names['ref_gapped'], names['query_gapped'])
    names['ref_unique'], names['query_unique'] = chain.print_only_unique(names['query_gapped'], names['ref_gapped'])

    current_file = names['query_gapped']
    num_lines = len(chain.alignment)
    make_image_from_repeat_fasta(current_file, consensus_width, num_lines)


def make_image_from_repeat_fasta(current_file, consensus_width, num_lines):
    layout = TransposonLayout()
    output_folder, output_file_name = os.path.split(current_file)
    layout.create_image_from_preprocessed_alignment(current_file, consensus_width, num_lines, output_folder, output_file_name)


if __name__ == '__main__':
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(BASE_DIR)
    # align_annotation('data\\RepeatMasker_all_alignment.csv',
    #                  'data\\hg38.fa',
    #                  'data\\panTro4_chr20.fa',
    #                  'data\\hg38ToPanTro4.over.chain')
    make_image_from_repeat_fasta('chr1\\panTro4_to_hg38_chr1_gapped.fa', 972, 107 * 2)