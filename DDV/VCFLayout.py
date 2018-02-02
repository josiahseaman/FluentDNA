from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import os
import traceback
from datetime import datetime

import math
from PIL import Image, ImageDraw
from natsort import natsorted

from DNASkittleUtils.CommandLineUtils import just_the_name
from DNASkittleUtils.Contigs import Contig, chunks

from DDV import gap_char
from DDV.DDVUtils import LayoutLevel
from DDV.TileLayout import TileLayout


def read_vcf_to_contig(vcf_file_path):
    """Equivalent to DNASkittleUtils.Contigs.read_contigs but for VCF SNP files.
    Reads the sequence of detected SNPS and appends them together as if they were
    a sequence.  Each specimen gets on 'seq' that is its unique set of SNPs."""
    import vcf
    samples = []
    vcf_reader = vcf.Reader(open(vcf_file_path, 'r'))
    for record in vcf_reader:
        if not samples:
            samples = [list() for i in range(2 * len(record.samples))]
        for i, specimen in enumerate(record.samples):
            bases = specimen.gt_bases
            bases = '..' if bases is None else bases.replace('/', '')
            samples[i * 2].append(bases[0])  # each phase of the diploid goes in a separate sequence
            samples[i * 2 + 1].append(bases[1])
    contigs = [Contig(name=vcf_reader.samples[i//2],
                      seq=''.join(samples[i])) for i in range(len(samples))]
    return contigs


class VCFLayout(TileLayout):
    def __init__(self, **kwargs):
        super(VCFLayout, self).__init__(**kwargs)
        self.line_width = 1000  # size of the display line before wrapping to the next alignment row
        self.column_height = 1000


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
            block = self.rearrange_alignment_into_megarows(specimens, vcf_file_path)
            block.width = self.line_width
            block.height = len(specimens)  # 2 rows per samples specimen diploid
            self.contigs.append(block)
        return self.contigs


    def adjust_layout_with_new_row_height(self, height):
        self.column_height = height # number of specimens this file

        # noinspection PyListCreation
        self.levels = [
            LayoutLevel("XInColumn", self.line_width, 1, 0),  # [0]
            LayoutLevel("LineInColumn", self.column_height, self.line_width * 1, 0)  # [1]
        ]
        self.levels.append(LayoutLevel("ColumnInRow", 1, padding=6, levels=self.levels))  # [2]
        # TODO: modulo could be self.image.width // self.line_width
        self.levels.append(LayoutLevel("RowInTile", 1000, levels=self.levels))  # [3]

    # def initialize_image_by_sequence_dimensions(self, consensus_width=None, num_lines=None):
    #     consensus_width = sum([x.consensus_width + 20 for x in self.contigs]) // len(self.contigs)
    #     min_width = max([x.consensus_width for x in self.contigs])
    #     consensus_width = max(consensus_width, min_width)
    #     self.layout_based_on_repeat_size(consensus_width)
    #     heights = [len(x.seq) // x.consensus_width for x in self.contigs]
    #     num_lines = sum(heights) + 20 * len(self.contigs)
    #     self.adjust_layout_with_new_row_height(heights)  # max
    #     min_height = self.column_height
    #     print("Average Width", consensus_width, "Genes", num_lines)
    #     self.image_length = consensus_width * num_lines
    #     basic_dimension = int(math.sqrt(self.image_length))
    #     width = max(min_width, basic_dimension)
    #     height = max(min_height, basic_dimension)
    #     self.image = Image.new('RGB', (width, height), "white")
    #     self.draw = ImageDraw.Draw(self.image)
    #     self.pixels = self.image.load()
    #     print("Image is ", self.image.width, "x", self.image.height)

    def rearrange_alignment_into_megarows(self, specimens, vcf_file_path):
        """
        :type specimens list(Contig)
        :return: Contig containing the whole file reformatted for self.line_width layout
        """
        display_lines = []
        for specimen in specimens:
            # each specimen generates a handful of chunks, each chunk is 1000bp long
            sample_seq_chunks = [''.join(line) for line in chunks(specimen.seq, self.line_width)]
            remainder = self.line_width - len(sample_seq_chunks[-1])
            sample_seq_chunks[-1] = sample_seq_chunks[-1] + (gap_char * remainder)
            display_lines.append(sample_seq_chunks)
        alignment_chunks = list(zip(*display_lines))  # transpose chunk lines to align samples
        chunk_seqs = [''.join(seq) for seq in alignment_chunks] # lines joined into blocks, blocks
        final_seq = ''.join(chunk_seqs)  # blocks concatenated for layout to separate
        return Contig(just_the_name(vcf_file_path), final_seq)


    def draw_nucleotides(self):
        self.draw_nucleotides_in_variable_row_height()


    def draw_nucleotides_in_variable_row_height(self):
        """Layout a whole set of different repeat types with different heights, but a fixed width (1,000).
          Wrapping to the next mega column is determined by hitting the bottom of the allocated image."""
        for contig in self.contigs:
            assert contig.height, "You must set the height in order to use this layout"
            self.adjust_layout_with_new_row_height(contig.height)

            contig_progress = 0
            seq_length = len(contig.seq)
            line_width = contig.width
            for cx in range(0, seq_length, line_width):
                x, y = self.position_on_screen(contig_progress)
                if x + contig.width + 10 >= self.image.width:
                    print("Ran into bottom of image at", contig.name)
                    return contig  # can't fit anything more
                # if y + contig.height >= self.image.height:
                #     contig_progress = 0  # reset to top of image, new column
                #     self.skip_to_next_mega_row()
                #     x, y = self.position_on_screen(contig_progress)
                # if y == self.origin[1]:  # first line in a column
                #     self.draw_repeat_title(contig, x, y) TODO

                remaining = min(line_width, seq_length - cx)
                if x + remaining >= self.image.width:
                    print("Ran into right side of image at", contig.name)
                    return contig  # can't fit anything more
                contig_progress += remaining
                for i in range(remaining):
                    nuc = contig.seq[cx + i]
                    # if nuc != gap_char:
                    self.draw_pixel(nuc, x + i, y)

            print("Drew", contig.name, "at", self.position_on_screen(contig_progress))
            columns_consumed = int(math.ceil(contig_progress / self.levels[2].chunk_size))
            self.origin[0] += columns_consumed * self.levels[2].thickness
        print('')


    def skip_to_next_mega_row(self):
        print("Skipping to next row:", self.origin)
        self.origin[0] = self.levels[2].padding  # start at left again
        self.origin[1] += self.levels[3].thickness  # go to next mega row


    def initialize_image_by_sequence_dimensions(self):
        min_width = self.line_width
        min_height = max([x.height for x in self.contigs])
        self.image_length = sum([len(x.seq) for x in self.contigs])
        basic = int(math.sqrt(self.image_length))
        width, height = max(min_width, basic), max(min_height, basic)
        self.image = Image.new('RGB', (width, height), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.pixels = self.image.load()
        print("Image is ", self.image.width, "x", self.image.height)

    def contig_json(self):
        return '[]'

    def levels_json(self):
        return '[]'
