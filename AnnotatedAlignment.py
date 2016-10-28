import os
from array import array
from Annotations import create_fasta_from_annotation
from ChainParser import ChainParser
from DDVUtils import Batch, pluck_contig, write_complete_fasta, first_word


class AnnotatedAlignment(ChainParser):
    def __init__(self, chain_name,
                 first_source, first_annotation,
                 second_source, second_annotation,
                 output_prefix,
                 trial_run=False, separate_translocations=False, squish_gaps=False,
                 show_translocations_only=False, aligned_only=False):
        super(AnnotatedAlignment, self).__init__(chain_name, first_source, second_source,
                                                 output_prefix, trial_run=trial_run,
                                                 separate_translocations=separate_translocations,
                                                 squish_gaps=squish_gaps,
                                                 show_translocations_only=show_translocations_only,
                                                 aligned_only=aligned_only)
        self.ref_annotation_source = first_annotation
        self.query_annotation_source = second_annotation


    def gap_annotation_metadata(self):
        """Modifies the annotation meta data file to compensate for coordinate shifts"""
        for pair in self.alignment:
            # TODO: output headers JSON again, but with indices compensated for gaps
            pass


    def switch_sequences(self, query_name, query_strand):
        # TODO: something
        return True


    def _parse_chromosome_in_chain(self, chromosome_name) -> Batch:
        names, ref_chr = self.setup_for_reference_chromosome(chromosome_name)
        master_chain = [chain for chain in self.chain_list if chain.tName == ref_chr and chain.qName == ref_chr and chain.qStrand == '+']  # TODO: remove this line
        self.create_alignment_from_relevant_chains(ref_chr, relevant_chains=master_chain)

        self.ref_sequence = pluck_contig(ref_chr, self.ref_source)  # only need the reference chromosome read, skip the others
        self.query_sequence = self.query_contigs[ref_chr]  # TODO: remove this line
        self.create_fasta_from_composite_alignment()
        names['query_gapped'], names['ref_gapped'] = self.write_gapped_fasta(names['ref'], names['query'])
        self.query_seq_gapped = array('u', '')
        self.ref_seq_gapped = array('u', '')
        # At this point we have created two gapped sequence fastas

        # Now create two annotation fastas so that we can gap them
        ref_annotation_fasta = os.path.join(self.output_folder, first_word(self.ref_source) + '_annotation_' + ref_chr + '.fa')
        create_fasta_from_annotation(self.ref_annotation_source, ref_chr, ref_annotation_fasta, len(self.ref_sequence))
        query_annotation_fasta = os.path.join(self.output_folder, first_word(self.query_source) + '_annotation_' + ref_chr + '.fa')
        create_fasta_from_annotation(self.query_annotation_source, ref_chr, query_annotation_fasta, len(self.query_contigs[ref_chr]))

        self.ref_sequence = pluck_contig(ref_chr, ref_annotation_fasta)
        self.query_sequence = pluck_contig(ref_chr, query_annotation_fasta)
        self.create_fasta_from_composite_alignment()
        self.markup_annotation_differences()
        # TODO: self.gap_annotation_metadata()
        names['r_anno_gap'], names['q_anno_gap'] = self.write_gapped_fasta(ref_annotation_fasta, query_annotation_fasta, False)

        # NOTE: Order of these appends DOES matter!
        self.output_fastas.append(names['r_anno_gap'])
        self.output_fastas.append(names['q_anno_gap'])
        self.output_fastas.append(names['ref_gapped'])
        self.output_fastas.append(names['query_gapped'])

        if True:  # self.trial_run:  # these files are never used in the viz
            del names['ref']
            del names['query']
        batch = Batch(chromosome_name, self.output_fastas, self.output_folder)
        self.output_folder = None  # clear the previous value
        return batch


    def markup_annotation_differences(self):
        # TODO: are there any headers or comments to scan past?
        shared_length = min(len(self.ref_seq_gapped), len(self.query_seq_gapped))
        human_unique, chimp_unique, shared_transcription = 0, 0, 0
        human_exons, chimp_exons, shared_exons = 0, 0, 0
        for i in range(shared_length):
            R = self.ref_seq_gapped[i]
            Q = self.query_seq_gapped[i]
            if R == 'G':  # Set of conditions where either exons or genes disagree
                if Q != 'G':
                    self.query_seq_gapped[i] = 'C'  # mark query deficiencies in blue
                    human_exons += 1
                else:
                    shared_exons += 1
            elif Q == 'G':
                self.ref_seq_gapped[i] = 'A'  # mark query deficiencies in red
                chimp_exons += 1

            if R == 'T':
                if Q != 'T':
                    human_unique += 1
                else:
                    shared_transcription += 1
            elif Q == 'T':
                chimp_unique += 1

        print("Differences in annotations are marked in red.")
        total_transcribed = human_unique + chimp_unique + shared_transcription
        total_exons = human_exons + chimp_exons + shared_exons
        print("Reference: %.2f percent transcriptional overlap" % ((human_unique + shared_transcription) / total_transcribed))
        print("Query:     %.2f percent transcriptional overlap" % ((chimp_unique + shared_transcription) / total_transcribed))
        print("Reference: %.2f percent exon bp overlap" % ((human_exons + shared_exons) / total_exons))
        print("Query:     %.2f percent exon bp overlap" % ((chimp_exons + shared_exons) / total_exons))
        with open(os.path.join(self.output_folder, 'stats.txt'), '+') as stats:
            values = 'human_unique chimp_unique shared_transcription human_exons chimp_exons shared_exons'.split()
            for val in values:
                stats.write('%s\t%i\n' % (val, locals()[val]))


if __name__ == '__main__':
    output_name = 'hg38_panTro4_annotated_'
    base_path = os.path.join('.', 'www-data', 'dnadata', output_name)

    chimp_annotation = r'HongKong\PanTro_refseq2.1.4_genes.gtf'
    human_anno = r'HongKong\Hg38_genes.gtf'
    # aligner = AnnotatedAlignment('hg38ToPanTro4.over.chain', 'panTro4.fa', 'hg38.fa', 'hg38_panTro4_annotated_')
    aligner = AnnotatedAlignment('hg38ToPanTro4.over.chain', 'hg38.fa', human_anno, 'panTro4.fa', chimp_annotation, base_path)
    aligner.parse_chain(['chr20'])
    # Force stop after gapped fasta file is written
