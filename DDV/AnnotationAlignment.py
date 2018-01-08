import os
from array import array

from ChainParser import ChainParser
from DNASkittleUtils.CommandLineUtils import just_the_name
from RepeatAnnotations import max_consensus_width, read_repeatmasker_csv, \
    filter_repeats_by_chromosome_and_family
from Span import alignment_chopping_index, AlignedSpans, Span
from TransposonLayout import TransposonLayout

from DDVUtils import make_output_dir_with_suffix


# class AnnotationAlignment(object):
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
    chromosome = 'chr19'
    all_repeat_entries = read_repeatmasker_csv(annotation_filename, 'genoName', chromosome)
    repeat_families = set(x.rep_family.replace('?', '#') for x in all_repeat_entries if x.rep_family)

    chain = ChainParser(chain_file, ref_fasta, query_fasta, '')
    # We want to do pieces of the contents of _parse_chromosome_in_chain
    names, ref_chr = chain.setup_for_reference_chromosome(chromosome)  # TODO: make this the parent folder
    chain.create_alignment_from_relevant_chains(chromosome)

    for repeat_type_name in repeat_families:
        ending = chromosome + '_' + repeat_type_name
        repeat_entries = filter_repeats_by_chromosome_and_family(all_repeat_entries, chromosome, repeat_type_name)
        repeat_entries.sort(key=lambda x: -len(x))
        try:
            consensus_width = max_consensus_width(repeat_entries)
        except ValueError:
            continue  # There's no usable entries left to use

        # modify chain.alignment to only contain annotated stretches
        chain.output_folder = make_output_dir_with_suffix(chain.output_prefix, ending)  # create a folder specifically for this repeat
        chain.query_seq_gapped = array('u', '')  # these need to be cleared so they don't accumulate the previous family
        chain.ref_seq_gapped = array('u', '')
        trimmed_alignment = create_aligned_annotation_fragments(chain.alignment, repeat_entries)
        chain.create_fasta_from_composite_alignment(previous_chr=None, alignment=trimmed_alignment)  # populates chain.query_seq_gapped and ref_seq_gapped
        # or if no fasta output is wanted: query_uniq_array, ref_uniq_array = chain.compute_unique_sequence()

        names['ref_gapped'], names['query_gapped'] = chain.write_gapped_fasta(names['ref'], names['query'], prepend_output_folder=True)
        names['ref_unique'], names['query_unique'] = chain.print_only_unique(names['query_gapped'], names['ref_gapped'])

        num_lines = len(trimmed_alignment) + 30 + 10  # 30 for title, 10 for origin
        for key in ['ref_gapped', 'query_gapped', 'ref_unique', 'query_unique']:
            make_image_from_repeat_fasta(names[key], consensus_width, num_lines)


def make_image_from_repeat_fasta(current_file, consensus_width, num_lines):
    layout = TransposonLayout()
    output_folder = os.path.dirname(current_file)
    output_file_name = just_the_name(current_file)
    layout.create_image_from_preprocessed_alignment(current_file, consensus_width, num_lines, output_folder, output_file_name)


if __name__ == '__main__':
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(BASE_DIR)
    align_annotation('data\\RepeatMasker_all_alignment.csv',
                     'data\\hg38.fa',
                     'data\\panTro4.fa',
                     'data\\hg38ToPanTro4.over.chain')
    # make_image_from_repeat_fasta('chr1\\panTro4_to_hg38_chr1_gapped.fa', 972, 107 * 2)