#!/usr/bin/env python
"""
DEPRECATION WARNING: This module has received only limited testing with one dataset. If you
want to replicate similar results with your own data, you'll likely need to
modify this python file.
This module converts two annotations files into two FASTA files to be used as
a visualization track alongside
two genomes and uses a whole genome alignment .chain file to align the set
of four "tracks" together.  This allows you to compare two annotations in the
context of their respective sequences.
"""
import os
from DNASkittleUtils.DDVUtils import editable_str

from DDV.ChainParser import ChainParser, scan_past_header, Batch
from DNASkittleUtils.Contigs import pluck_contig
from DNASkittleUtils.DDVUtils import first_word, ReverseComplement

from DDV.Annotations import create_fasta_from_annotation, GFF


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
        self.annotation_phase = False
        self.query_GFF = GFF(self.query_annotation_source)



    def gap_annotation_metadata(self):
        """Modifies the annotation meta data file to compensate for coordinate shifts"""
        for pair in self.alignment:
            # TODO: output headers JSON again, but with indices compensated for gaps
            pass


    # def missing_query_sequence(self, query_name):
    #     annotation_fasta = create_fasta_from_annotation(self.query_GFF, query_name, None)
    #     if annotation_fasta:  # content
    #         self.query_contigs[query_name] = annotation_fasta
    #         self.query_sequence = annotation_fasta
    #         return True
    #     return super(AnnotatedAlignment, self).missing_query_sequence(query_name)


    def rev_comp_contig(self, query_name):
        return ReverseComplement(self.query_contigs[query_name], annotation=self.annotation_phase)


    def _parse_chromosome_in_chain(self, chromosome_name):# -> Batch:
        print("=== Begin Annotated Alignment ===")
        names, ref_chr = self.setup_for_reference_chromosome(chromosome_name)
        self.create_alignment_from_relevant_chains(ref_chr)

        self.ref_sequence = pluck_contig(ref_chr, self.ref_source)  # only need the reference chromosome read, skip the others
        self.query_sequence = self.query_contigs[ref_chr]  # TODO: remove this line
        self.create_fasta_from_composite_alignment()
        names['ref_gapped'], names['query_gapped'] = self.write_gapped_fasta(names['ref'], names['query'])

        self.query_seq_gapped = editable_str('')
        self.ref_seq_gapped = editable_str('')
        self.query_contigs = {}
        self.stored_rev_comps = {}
        self.annotation_phase = True
        # At this point we have created two gapped sequence fastas

        query_annotation_fasta, ref_annotation_fasta = self.load_annotation_fastas(ref_chr)

        self.create_fasta_from_composite_alignment(previous_chr=(ref_chr, '+'))
        self.markup_annotation_differences()
        # TODO: self.gap_annotation_metadata()
        names['r_anno_gap'], names['q_anno_gap'] = self.write_gapped_fasta(ref_annotation_fasta, query_annotation_fasta, False)
        self.write_stats_file()
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


    def load_annotation_fastas(self, ref_chr):
        # Now create two annotation fastas so that we can gap them
        #TODO: check for         copy_to_sources(self.output_folder, )
        self.ref_sequence = create_fasta_from_annotation(self.ref_annotation_source, [ref_chr], None)[0].seq
        self.query_contigs[ref_chr] = create_fasta_from_annotation(self.query_GFF, [ref_chr], None)[0].seq

        # TODO: these four lines are place holders for tracking real annotation ref and query contigs
        # self.ref_sequence = pluck_contig(ref_chr, ref_annotation_fasta)
        # self.query_sequence = pluck_contig(ref_chr, query_annotation_fasta)
        self.query_sequence = self.query_contigs[ref_chr]

        ref_annotation_fasta = os.path.join(self.output_folder, first_word(self.ref_source) + '_annotation_' + ref_chr + '.fa')
        query_annotation_fasta = os.path.join(self.output_folder, first_word(self.query_source) + '_annotation_' + ref_chr + '.fa')
        return query_annotation_fasta, ref_annotation_fasta


    def markup_annotation_differences(self):
        print("Marking annotation differences...")
        r = scan_past_header(self.ref_seq_gapped, 0)
        q = scan_past_header(self.query_seq_gapped, 0)
        while q < len(self.query_seq_gapped) and r < len(self.ref_seq_gapped):
            # query_uniq_array is already initialized to contain header characters if separating translocations
            r = scan_past_header(self.ref_seq_gapped, r, )  # take_shortcuts=not self.separate_translocations)
            q = scan_past_header(self.query_seq_gapped, q, )  # take_shortcuts=not self.separate_translocations)
            R = self.ref_seq_gapped[r]
            Q = self.query_seq_gapped[q]
            if R == 'G':  # Set of conditions where either exons or genes disagree
                if Q != 'G':
                    self.query_seq_gapped[q] = 'C'  # mark query deficiencies in blue
                    self.stats['human_exons_bp'] += 1
                else:
                    self.stats['shared_exons_bp'] += 1
            elif Q == 'G':
                self.ref_seq_gapped[r] = 'A'  # mark query deficiencies in red
                self.stats['chimp_exons_bp'] += 1

            if R == 'T':
                if Q != 'T':
                    self.stats['human_unique_transcription_bp'] += 1
                else:
                    self.stats['shared_transcription_bp'] += 1
            elif Q == 'T':
                self.stats['chimp_unique_transcription_bp'] += 1
            q += 1
            r += 1

        print("Differences in annotations are marked in red.")
        self.stats['total_bp_transcribed'] = self.stats['human_unique_transcription_bp'] + self.stats['chimp_unique_transcription_bp'] + self.stats['shared_transcription_bp']
        self.stats['total_exons_bp'] = self.stats['human_exons_bp'] + self.stats['chimp_exons_bp'] + self.stats['shared_exons_bp']
        print("Of all transcription area: %.1f%% is uniquely chimp, %.1f%% is uniquely human, %.1f%% is shared" %
              (self.stats['chimp_unique_transcription_bp'] / self.stats['total_bp_transcribed'] * 100,
               self.stats['human_unique_transcription_bp'] / self.stats['total_bp_transcribed'] * 100,
               self.stats['shared_transcription_bp'] / self.stats['total_bp_transcribed'] * 100))
        print("Of all Exons bp: %.1f%% is uniquely chimp, %.1f%% is uniquely human %.1f%% is shared" %
              (self.stats['chimp_exons_bp'] / self.stats['total_exons_bp'] * 100,
               self.stats['human_exons_bp'] / self.stats['total_exons_bp'] * 100,
               self.stats['shared_exons_bp'] / self.stats['total_exons_bp'] * 100))




if __name__ == '__main__':
    output_name = 'hg38_panTro4_annotated_'
    base_path = os.path.join('.', 'www-data', 'dnadata', output_name)
    chimp_annotation = r'data\PanTro_refseq2.1.4_genes.gtf'
    human_anno = r'data\Hg38_genes.gtf'
    aligner = AnnotatedAlignment('hg38ToPanTro4.over.chain', 'hg38.fa', human_anno, 'panTro4.fa', chimp_annotation, base_path)
    aligner.parse_chain(['chr20'])

    #### ==== Command Line Configuration === ####
    # fluentdna.py --chainfile=hg38ToPanTro4.over.chain --fasta=hg38.fa --extrafastas panTro4.fa --ref_annotation=DDV\\data\Hg38_genes.gtf
    # --query_annotation=DDV\data\PanTro_refseq2.1.4_genes.gtf --outname=hg38_panTro4_annotated_
