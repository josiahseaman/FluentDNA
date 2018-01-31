from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import os
import traceback
from datetime import datetime

from DNASkittleUtils.CommandLineUtils import just_the_name
from DNASkittleUtils.Contigs import Contig
from natsort import natsorted

from DDV.MultipleAlignmentLayout import MultipleAlignmentLayout


def read_vcf_to_contig(vcf_file_path):
    """Equivalent to DNASkittleUtils.Contigs.read_contigs but for VCF SNP files.
    Reads the sequence of detected SNPS and appends them together as if they were
    a sequence.  Each specimen gets on 'seq' that is its unique set of SNPs."""
    import vcf
    samples = []
    vcf_reader = vcf.Reader(open(vcf_file_path, 'r'))
    for record in vcf_reader:
        if not samples:
            samples = [list() for i in record.samples]
        for i, specimen in enumerate(record.samples):
            bases = specimen.gt_bases
            bases = '..' if bases is None else bases.replace('/', '')
            samples[i].append(bases)
    contigs = [Contig(name=vcf_reader.samples[i],
                      seq=''.join(samples[i])) for i in range(len(samples))]
    return contigs


class VCFLayout(MultipleAlignmentLayout):
    def __init__(self):
        super(VCFLayout, self).__init__()


    def process_all_alignments(self, input_vcf_folder, output_folder, output_file_name):
        self.origin[1] += self.levels[5].padding  # One full Row of padding for Title
        start_time = datetime.now()
        self.translate_vcf_folder_to_contigs(input_vcf_folder)
        print("Converted VCF Files :", datetime.now() - start_time)

        self.initialize_image_by_sequence_dimensions()
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)


    def translate_vcf_folder_to_contigs(self, input_vcf_folder):
        from glob import glob
        self.contigs = []
        for vcf_file_path in natsorted(glob(os.path.join(input_vcf_folder, '*.vcf'))):
            specimens = read_vcf_to_contig(vcf_file_path)
            consensus_width = max([len(x.seq) for x in specimens])
            block = Contig(just_the_name(vcf_file_path),
                           ''.join([x.seq for x in specimens]),)
            block.consensus_width = consensus_width
            self.contigs.append(block)
        return self.contigs

    def set_column_height(self, heights):
        self.column_height = max(heights)  # For clusters of VCF, there's likely to be the same
        # number of specimens for every file
