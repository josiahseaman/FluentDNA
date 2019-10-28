from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import os
import shutil
import traceback
from datetime import datetime

from DNASkittleUtils.Contigs import read_contigs
from PIL import Image, ImageDraw

import math
from FluentDNA.TileLayout import hex_to_rgb, TileLayout, is_protein_sequence
from natsort import natsorted

from FluentDNA.FluentDNAUtils import make_output_directory
from FluentDNA.Layouts import level_layout_factory


def fastas_in_folder(input_fasta_folder):
    from glob import glob
    # If I was actually given a file name and not a directory, just process the one file
    if not os.path.isdir(input_fasta_folder) and os.path.exists(input_fasta_folder):
        return [input_fasta_folder]
    else:
        return list(natsorted(glob(os.path.join(input_fasta_folder, '*.fa*'))))



class MultipleAlignmentLayout(TileLayout):
    def __init__(self, sort_contigs=False, **kwargs):
        kwargs['low_contrast'] = True
        kwargs['sort_contigs'] = True
        super(MultipleAlignmentLayout, self).__init__(**kwargs)
        self.all_contents = {}  # (filename: contigs) output_fasta() determines order of fasta_sources
        self.current_column_height = 20
        self.next_origin = [self.border_width, 30] # margin for titles, incremented each MSA
        self.single_file = False  # flag to detect a single, long MSA
        self.sort_contigs = sort_contigs
        self.title_height_px = 10
        self.x_pad = 20  # whitespace between MSA blocks
        self.y_pad = 20  # vertical pad to include space for titles

        #### Rasmol 'Amino' Protein colors
        # self.palette['A'] = hex_to_rgb('C8C8C8')
        # self.palette['R'] = hex_to_rgb('145AFF')
        # self.palette['N'] = hex_to_rgb('00DCDC')
        # self.palette['D'] = hex_to_rgb('E60A0A')
        # self.palette['B'] = hex_to_rgb('E6E600')
        # self.palette['C'] = hex_to_rgb('00DCDC')
        # self.palette['E'] = hex_to_rgb('E60A0A')
        # self.palette['Q'] = hex_to_rgb('EBEBEB')
        # self.palette['Z'] = hex_to_rgb('8282D2')
        # self.palette['G'] = hex_to_rgb('0F820F')
        # self.palette['H'] = hex_to_rgb('0F820F')
        # self.palette['I'] = hex_to_rgb('145AFF')
        # self.palette['L'] = hex_to_rgb('E6E600')
        # self.palette['K'] = hex_to_rgb('3232AA')
        # self.palette['M'] = hex_to_rgb('DC9682')
        # self.palette['F'] = hex_to_rgb('FA9600')
        # self.palette['P'] = hex_to_rgb('FA9600')
        # self.palette['S'] = hex_to_rgb('B45AB4')
        # self.palette['T'] = hex_to_rgb('3232AA')
        # self.palette['W'] = hex_to_rgb('0F820F')
        # self.palette['Y'] = hex_to_rgb('FF69B4')
        # self.palette['V'] = hex_to_rgb('FF69B4')
        # self.palette['X'] = hex_to_rgb('FF6100')


        #### Rasmol Shapely Protein colors
        # Assert full Rasmol colors over Nucleotide values (different from TileLayout)
        # self.palette['A'] = hex_to_rgb('8CFF8C')
        # self.palette['R'] = hex_to_rgb('00007C')
        # self.palette['N'] = hex_to_rgb('FF7C70')
        # self.palette['D'] = hex_to_rgb('A00042')
        # self.palette['B'] = hex_to_rgb('FFFF70')
        # self.palette['C'] = hex_to_rgb('FF4C4C')
        # self.palette['E'] = hex_to_rgb('660000')
        # self.palette['Q'] = hex_to_rgb('FFFFFF')
        # self.palette['Z'] = hex_to_rgb('7070FF')
        # self.palette['G'] = hex_to_rgb('004C00')
        # self.palette['H'] = hex_to_rgb('455E45')
        # self.palette['I'] = hex_to_rgb('4747B8')
        # self.palette['L'] = hex_to_rgb('B8A042')
        # self.palette['K'] = hex_to_rgb('534C52')
        # self.palette['M'] = hex_to_rgb('525252')
        # self.palette['F'] = hex_to_rgb('FF7042')
        # self.palette['P'] = hex_to_rgb('B84C00')
        # self.palette['S'] = hex_to_rgb('4F4600')
        # self.palette['T'] = hex_to_rgb('8C704C')
        # self.palette['W'] = hex_to_rgb('FF8CFF')
        # self.palette['Y'] = hex_to_rgb('FF00FF')
        # self.palette['V'] = hex_to_rgb('FF00FF')


    def process_all_alignments(self, input_fasta_folder, output_folder, output_file_name):
        start_time = datetime.now()
        make_output_directory(output_folder)
        self.preview_all_files(input_fasta_folder)
        self.calculate_mixed_layout()
        print("Tallied all contigs :", datetime.now() - start_time)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        #TODO: sort all layouts with corresponding sequence?

        for file_no, single_MSA in enumerate(self.fasta_sources):
            self.i_layout = file_no
            self.contigs = self.all_contents[single_MSA]
            # self.read_contigs_and_calc_padding(single_MSA, None)
            try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
                self.draw_nucleotides()
                if self.use_titles and not self.single_file:
                    self.draw_titles()
            except Exception as e:
                print('Encountered exception while drawing nucleotides:', '\n')
                traceback.print_exc()
            input_path = os.path.join(input_fasta_folder, single_MSA)
            self.output_fasta(output_folder, input_path, False, None, False,
                              append_fasta_sources=False, create_source_download=False)
        print("\nDrew Nucleotides:", datetime.now() - start_time)
        self.output_image(output_folder, output_file_name, False)
        print("Output Image in:", datetime.now() - start_time)
        target_folder = os.path.join(output_folder, 'sources', os.path.basename(input_fasta_folder))
        if not os.path.exists(target_folder):
            print("Copying entire sources directory:", input_fasta_folder)
            shutil.copytree(input_fasta_folder,
                            target_folder,
                            ignore=lambda src, names: [n for n in names if '.fa' not in n],
                            symlinks=True, )


    def draw_nucleotides(self, verbose=False):
        """Layout a whole set of different repeat types with different widths.  Column height is fixed,
        but column width varies constantly.  Wrapping to the next row is determined by hitting the
        edge of the allocated image."""
        super(MultipleAlignmentLayout, self).draw_nucleotides(verbose)


    def calc_all_padding(self):
        total_progress = 0
        seq_start, title_length = 0, 0
        widest_sequence = 0
        for i, contig in enumerate(self.contigs):  # Type: class DNASkittleUtils.Contigs.Contig
            contig.reset_padding = 0
            contig.tail_padding = 0
            widest_sequence = max(widest_sequence, len(contig.seq))
            contig.consensus_width = widest_sequence
            #First contig of each MSA has a 10px tall title
            contig.title_padding = 0
            if i == 0 and not self.single_file:
                contig.title_padding = widest_sequence * self.title_height_px
            contig.nuc_title_start = seq_start
            contig.nuc_seq_start = seq_start + title_length
            #at the moment these values are the same but they have different meanings
            total_progress += len(contig.seq) + contig.title_padding  # pointer in image
            seq_start += title_length + len(contig.seq)  # pointer in text
        return total_progress



    def prepare_image(self, image_length, width=None, height=None):
        if not width or not height:
            width, height = self.max_dimensions(image_length)
        print("Image dimensions are", width, "x", height, "pixels")
        self.image = Image.new(self.pil_mode, (width, height), hex_to_rgb('#FFFFFF'))#ui_grey)
        self.draw = ImageDraw.Draw(self.image)
        self.pixels = self.image.load()

    def guess_image_dimensions(self):
        margin = self.border_width*2
        max_w = max([x.consensus_width for source in self.all_contents.values() for x in source]) + margin
        max_h = max([len(source) + self.title_height_px for source in self.all_contents.values()]) + margin
        #TODO check if max_h is too large and needs to be interrupted by layout: one full row
        areas = []
        for source in self.all_contents.values():
            areas.append((source[-1].consensus_width + self.x_pad) * (len(source) + self.y_pad))
        area = sum(areas) if not self.single_file else len(self.contigs[0].seq) * len(self.all_contents) *1.2
        self.image_length = int(area * 1.2)
        square_dim = int(math.sqrt(self.image_length))
        desired_width = 5 * square_dim // 3
        # TODO still not a great algorithm
        image_wh = [desired_width, max(max_h, self.image_length // desired_width)]
        return image_wh

    def calculate_mixed_layout(self):
        """Do complete layout of image, then decide its dimensions
        All the layout brains go here"""
        #TODO: sort contigs
        # bad_contigs = [c for c in self.contigs if not c.consensus_width]
        # for contig in bad_contigs:
        #     print("Error while reading FASTA. Skipping: %s" % contig.name)
        #     self.contigs.remove(contig)
        # if self.sort_contigs:
        #     self.contigs.sort(key=lambda x: -x.height)

        image_wh = self.guess_image_dimensions()

        #unsorted, largest height per row, tends to be less dense
        self.each_layout = []  # delete old default layout
        for filename in self.fasta_sources:
            source = self.all_contents[filename]
            height = len(source) + self.title_height_px
            width = source[0].consensus_width
            self.layout_based_on_repeat_size(width, height, image_wh[0], source)
        self.i_layout = 0 # drawing starts at the beginning

        adjusted_height = self.next_origin[1] + self.current_column_height + self.y_pad  #could extend image
        if self.single_file:
            adjusted_height = image_wh[1]
        self.prepare_image(0, image_wh[0], adjusted_height)



    def preview_all_files(self, input_fasta_folder):
        """Populates fasta_sources with files from a directory"""
        files = fastas_in_folder(input_fasta_folder)
        self.single_file = len(files) == 1
        if self.single_file:
            self.spread_large_MSA_source(files[0])
        else:
            for single_MSA in files:
                self.read_contigs_and_calc_padding(single_MSA, None)
                fasta_name = os.path.basename(single_MSA)
                self.fasta_sources.append(fasta_name)
                self.all_contents[fasta_name] = self.contigs  # store contigs so the can be wiped
            if self.sort_contigs:  # do this before self.each_layout is created in order
                heights = [(len(self.all_contents[fasta_name]), fasta_name) for fasta_name in self.fasta_sources]
                heights.sort(key=lambda pair: -pair[0])  # largest number of sequences first
                self.fasta_sources = [pair[1] for pair in heights]  # override old ordering


    def layout_based_on_repeat_size(self, width, height, max_width, contigs):
        """change layout to match dimensions of the repeat
        """
        total_width = width + self.x_pad
        max_rows = 1000
        if height > max_rows :
            columns = math.ceil(height / max_rows)
            total_width = (width + self.x_pad) * columns
            height = min(max_rows, height)
        if self.single_file:
            self.layout_phased_file(width, height, max_width)
        else:  # Typical case with many small MSA
            # skip to next mega row
            if self.next_origin[0] + total_width - self.x_pad + 1 >= max_width:
                self.next_origin[0] = self.border_width
                self.next_origin[1] += self.current_column_height + self.y_pad
                self.current_column_height = 1  # reset

            self.current_column_height = max(height, self.current_column_height)
            modulos = [width, height, 9999, 9999]
            padding = [0, 0, self.x_pad, self.x_pad * 3]
            self.each_layout.append(level_layout_factory(modulos, padding, self.next_origin))
            self.next_origin[0] += total_width  # scoot next_origin by width we just used up
        self.i_layout = len(self.each_layout) - 1  # select current layout


    def layout_phased_file(self, width, height, max_width):
        self.each_layout = []
        usable_width = min(width, max_width - (self.border_width * 2))
        # TODO more than one large MSA
        height = len(self.fasta_sources)  # number of individuals
        padding_between_mega_rows = 1
        n_rows = math.ceil(len(self.contigs[0].seq) / usable_width)
        modulos = [usable_width, 1, 1, n_rows]
        padding = [0, height + padding_between_mega_rows, 0, height + padding_between_mega_rows]
        for y, row in enumerate(self.fasta_sources):  # one layout for mouse over of each individual
            self.next_origin[0] = self.border_width
            self.next_origin[1] = 30 + y
            self.each_layout.append(level_layout_factory(modulos, padding, self.next_origin))
        # move origin to bottom of image
        self.next_origin[1] += n_rows * \
                               (height + padding_between_mega_rows)
        self.i_layout = len(self.each_layout) - 1  # select current layout


    def draw_titles(self):
        """Draw one title for each file (MSA Block) that includes many contigs.
        We use a fake contig with no sequence and the name of file, instead of the name of the first
        line of the MSA."""
        filename = os.path.basename(self.fasta_sources[self.i_layout])
        contig_name = os.path.splitext(filename)[0]  # should be basename
        upper_left = list(self.position_on_screen(0))
        upper_left[1] -= 2  # bit of margin, should be less than self.border_width
        font_size = 9  # font sizes: [9, 38, 380, 380 * 2]
        title_width = self.levels.base_width // 6
        title_lines = 1
        self.write_title(contig_name, self.levels.base_width, self.title_height_px, font_size,
                         title_lines, title_width, upper_left, False, self.image)

    def spread_large_MSA_source(self, fasta_path):
        individuals = read_contigs(fasta_path)
        self.contigs = individuals
        self.fasta_sources = [os.path.basename(fasta_path) + str(i) for i in range(len(individuals))]
        self.all_contents = {source: [individuals[i]] for i, source in enumerate(self.fasta_sources)}
        self.protein_palette = is_protein_sequence(self.contigs[0])

        # Zero padding
        for name, container in self.all_contents.items():
            contig = container[0]
            contig.reset_padding = 0
            contig.title_padding = 0
            contig.tail_padding = 0
            contig.nuc_title_start = 0
            contig.nuc_seq_start = 0
            contig.consensus_width = len(contig.seq)
