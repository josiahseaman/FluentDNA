from ChainParser import ChainParser
from Span import alignment_chopping_index, AlignedSpans, Span
from DDVUtils import Contig
from RepeatAnnotations import max_consensus_width, read_repeatmasker_csv
from TransposonLayout import TransposonLayout, grab_aligned_repeat, filter_repeats_by_chromosome


# class AnnotationAlignment:
#     def __init__(self):


def remove_all_unannotated_sections(alignment, repeat_entries):
    aligned_repeats = []  # list of AlignedSpans
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
                # TODO: add gaps based on repStart and repEnd
                aligned_repeats.append(AlignedSpans(ref, query, 0, 0, current.is_master_chain))
            else:
                discarded += 1
    print("Discarded %i out of %i partial overlaps" % (discarded, len(repeat_entries)))
    return aligned_repeats


def align_annotation(annotation_filename, ref_fasta, query_fasta, chain_file):
    # We want to do pieces of the contents of _parse_chromosome_in_chain
    chain = ChainParser(chain_file, ref_fasta, query_fasta, 'data/repeat_test/',
                        aligned_only=True)  # This allows us to skip over all intervening sequences
    names, ref_chr = chain.setup_for_reference_chromosome('chr1')
    chain.create_alignment_from_relevant_chains(ref_chr)

    # modify chain.alignment to only contain annotated stretches
    repeat_entries = read_repeatmasker_csv(annotation_filename, 'repName', 'L1PA3')
    chain.alignment = remove_all_unannotated_sections(chain.alignment, repeat_entries)
    chain.create_fasta_from_composite_alignment()

    names['ref_gapped'], names['query_gapped'] = chain.write_gapped_fasta(names['ref'], names['query'])
    names['ref_unique'], names['query_unique'] = chain.print_only_unique(names['query_gapped'], names['ref_gapped'])
    # padding with 'X' in fasta

    #
    # layout = TransposonLayout()
    # layout.read_all_files(ref_fasta, annotation_filename, column='repName', rep_name='L1PA3')
    # processed_contigs = create_repeat_fasta_contigs(layout.repeat_entries, layout.contigs)
    print("Consensus Width: ", max_consensus_width(repeat_entries))


def create_repeat_fasta_contigs(repeat_entries, contigs):
    consensus_width = max_consensus_width(repeat_entries)
    processed_contigs = []
    for contig in contigs:
        ref_samples = []
        query_samples = []
        reps_on_chr = filter_repeats_by_chromosome(repeat_entries, contig.name)
        for fragment in reps_on_chr:
            ref_line = grab_aligned_repeat(consensus_width, contig, fragment)
            ref_samples.append(''.join(ref_line))

            query_fragment = fragment  # TODO: with alignment
            query_line = grab_aligned_repeat(consensus_width, contig, query_fragment)
            query_samples.append(''.join(query_line))

        processed_seq = ''.join(ref_samples)
        processed_contigs.append(Contig(contig.name, processed_seq, 0, 0, 0,
                                        0, 0))  # TODO: title_length currently doesn't have a title and might break mouse tracking
        # TODO: fasta for query as well
    return processed_contigs


