import os

from Annotations import create_fasta_from_annotation
from ChainParser import ChainParser
from DDVUtils import Batch, pluck_contig, write_complete_fasta


class AnnotatedAlignment(ChainParser):
    def __init__(self, chain_name, first_source, second_source, output_prefix,
                 trial_run=False, separate_translocations=False, squish_gaps=False,
                 show_translocations_only=False, aligned_only=False):
        super(AnnotatedAlignment, self).__init__(chain_name, first_source, second_source,
                                                 output_prefix, trial_run=trial_run,
                                                 separate_translocations=separate_translocations,
                                                 squish_gaps=squish_gaps,
                                                 show_translocations_only=show_translocations_only,
                                                 aligned_only=aligned_only)


    def create_fasta_from_composite_alignment(self):
        """Modifies the annotation meta data file to compensate for coordinate shifts"""
        super(AnnotatedAlignment, self).create_fasta_from_composite_alignment()
        for pair in self.alignment:
            # TODO: output headers JSON again, but with indices compensated for gaps
            pass


    def write_gapped_annotation(self, reference, ref_gap_annotation, query, query_gap_annotation):
        query_gap_name, ref_gap_name = self.write_gapped_fasta(reference, query)
        write_complete_fasta(ref_gap_annotation, self.ref_anno_gapped)
        write_complete_fasta(query_gap_annotation, self.query_anno_gapped)
        return query_gap_name, ref_gap_name


    def _parse_chromosome_in_chain(self, chromosome_name) -> Batch:
        names, ref_chr = self.setup_for_reference_chromosome(chromosome_name)

        self.ref_sequence = pluck_contig(ref_chr, self.ref_source)  # only need the reference chromosome read, skip the others

        self.create_alignment_from_relevant_chains()
        self.create_fasta_from_composite_alignment()
        ref_gap_annotation = 'hg38_gapped_annotation_chr20.fa'
        query_gap_annotation = 'panTro4_gapped_annotation_chr20.fa'
        names['query_gapped'], names['ref_gapped'] = self.write_gapped_annotation(names['ref'], ref_gap_annotation, names['query'], query_gap_annotation)
        # NOTE: Order of these appends DOES matter!
        self.output_fastas.append(ref_gap_annotation)
        self.output_fastas.append(query_gap_annotation)
        self.output_fastas.append(names['ref_gapped'])
        self.output_fastas.append(names['query_gapped'])
        print("Finished creating gapped fasta files", names['ref'], names['query'])

        if True:  # self.trial_run:  # these files are never used in the viz
            del names['ref']
            del names['query']
        batch = Batch(chromosome_name, self.output_fastas, self.output_folder)
        self.output_folder = None  # clear the previous value
        return batch


if __name__ == '__main__':
    annotation_file = r'HongKong\Pan_Troglodytes_refseq2.1.4.gtf'
    target_chromosome = 'chr20'
    chimp_anno_fasta = 'panTro4_annotation_' + target_chromosome + '.fa'
    create_fasta_from_annotation(annotation_file, target_chromosome, 63 * 1000 * 1000, chimp_anno_fasta)

    # aligner = AnnotatedAlignment('hg38ToPanTro4.over.chain', 'panTro4.fa', 'hg38.fa', 'hg38_panTro4_annotated_')
    aligner = AnnotatedAlignment('hg38ToPanTro4.over.chain', chimp_anno_fasta, 'hg38.fa', 'hg38_panTro4_annotated_')
    aligner.parse_chain(['chr20'])
    # Force stop after gapped fasta file is written
