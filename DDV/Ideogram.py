#!/usr/bin/env python
"""Plot a DNA sequence along a murray space filling curve.

This FluentDNA Python class is based on Yan Wong's code (unpublished) which is in turn
based on: algol version from https://alb.host.cs.st-andrews.ac.uk/cole/code.html
(see also http://www.aleph.se/andart/archives/2013/10/murray_polygons.html)

The murray polygon is stretched by 3 times in the X and 5 times in the Y direction
to improve the locality properties (ideally this should be 1:sqrt(3) = 0.577 rather
than 3:5 = 0.6) see https://arxiv.org/pdf/0806.4787.pdf

try e.g. python3 Ideogram.py -x 3 3 3 -y 3 3 3
"""
import sys
import traceback
from datetime import datetime
from itertools import chain

from DNASkittleUtils.Contigs import read_contigs

from DDV.DDVUtils import beep, make_output_directory
from DDV.HighlightedAnnotation import HighlightedAnnotation
import os
import numpy as np
from functools import reduce

from Layouts import LayoutFrame, LayoutLevel
from MultipleAlignmentLayout import fastas_in_folder
from TileLayout import hex_to_rgb


class IdeogramCoordinateFrame(LayoutFrame):
    def __init__(self, x_radices, y_radices, x_scale, y_scale, border_width):
        self.point_mapping = [] # for annotation and testing purposes
        self.fibre_padding = 3
        self.x_radices = x_radices
        self.y_radices = y_radices
        self.x_scale = x_scale
        self.y_scale = y_scale

        #levels creation
        modulos = list(chain(*zip(x_radices, y_radices)))
        padding = [0] * len(modulos)
        padding[-1] = self.fibre_padding  # very last y layer has padding to match hacked_padding()
        # This has the side effect that coordinates are wrong for the last half of the last
        # block in each row where hacked padding comes in the middle and mouse padding is at the end
        levels = [LayoutLevel(modulos[0], 1, padding[0]),
                  LayoutLevel(modulos[1], modulos[0], padding[1])]
        for i in range(2, len(modulos)):
            levels.append(LayoutLevel(modulos[i], padding=padding[i], levels=levels))
        self.levels = levels

        super(IdeogramCoordinateFrame, self).__init__((border_width, border_width), self.levels)  # levels is our iterable
        #Itering levels is less helpful in Ideogram than it is in regular LayoutFrame


    @property
    def base_width(self):
        """Shorthand for the column width value that is used often.  This can change
        based on the current self.i_layout."""
        return reduce(int.__mul__, self.x_radices)


    def build_coordinate_mapping(self, sequence_length):
        while self.levels[-1].chunk_size > sequence_length:
            self.levels = self.levels[:-1]  # drop unnecessary levels used in mouse calculation

        ndim_x = len(self.x_radices)
        ndim_y = len(self.y_radices)
        max_dim = max(ndim_x, ndim_y)

        radices = np.ones((max_dim, 2), dtype=np.int)
        digits = np.zeros((max_dim * 2), dtype=np.int)
        parities = np.ones((max_dim + 1, 2), dtype=np.int)
        curr_pos = np.array((0,0), dtype=np.int)  # stores y,x
        radices[0:ndim_x, 0] = self.x_radices
        radices[0:ndim_y, 1] = self.y_radices

        no_pts = np.prod(radices)
        radices.shape = np.prod(radices.shape)  # flatten

        if self.x_scale == 1 and self.y_scale == 1:
            self.build_loop_optimized(curr_pos, digits, no_pts, parities, radices, sequence_length)
        else:
            prevprev_pos = np.zeros((2,), dtype=np.int) # must start at 0 because of scale multiplaction
            prev_pos = np.zeros((2,), dtype=np.int)  # must start at 0 because of scale multiplaction
            curr_pos = np.zeros((2,), dtype=np.int)
            raise NotImplementedError("Scales beyond 1,1 are not currently implemented")
            # self.draw_loop_any_scale(curr_pos, digits, no_pts, parities, points_file,
            #                          prev_pos, prevprev_pos, radices, seq_iter, self.x_scale, self.y_scale)


    def build_loop_optimized(self, curr_pos, digits, no_pts, parities, radices, sequence_length):
        max_x = reduce(int.__mul__, self.x_radices) - 1 #+ self.origin[0]
        min_x = 0  #self.origin[0]
        odd = 0
        for pts in range(min(sequence_length, no_pts - 1)):
            place = increment(digits, radices, 0)
            parities[0:(place // 2 + 1), place % 2] *= -1
            place += 1
            odd = self.hacked_padding(curr_pos, min_x, max_x, odd, place)
            x = int(curr_pos[1])
            y = int(curr_pos[0])
            curr_pos[place % 2] += parities[place // 2, place % 2]
            self.point_mapping.append((x,y))


    def hacked_padding(self, curr_pos, min_x, max_x, odd, place):
        if place % 2 == 0:  # this is an y increments
            if place // 2 == len(self.x_radices) - 1:
                if curr_pos[1] == max_x or curr_pos[1] == min_x:
                    if odd == 1:
                        curr_pos[0] += self.fibre_padding  # y coordinates are in [0]
                    odd = (odd + 1) % 2
        return odd


    def position_on_screen(self, progress):
        """WARNING: This will not work until after self.draw_loop_optimized
         has populated self.point_mapping"""
        x, y = self.point_mapping[progress]
        return x + self.origin[0], y + self.origin[1]

    def relative_position(self, progress):
        return self.point_mapping[progress]


    def handle_multi_column_annotations(self, start, stop):
        """In 2D fractal layout, this method is much simpler since there's no columns per se.
        It may be a good idea to identify the largest chromatin fibre and ensure that labels
        don't straddle that boundary."""
        pts = self.point_mapping[start:stop]  # all the coordinates annotated by this region
        left, right = min(pts, key=lambda p: p[0])[0], max(pts, key=lambda p: p[0])[0]
        top, bottom = min(pts, key=lambda p: p[1])[1], max(pts, key=lambda p: p[1])[1]
        height = bottom - top
        width = right - left
        return width, height, left + self.origin[0], right + self.origin[0],\
               top + self.origin[1], bottom + self.origin[1]




class Ideogram(HighlightedAnnotation):
    def __init__(self, radix_settings, ref_annotation=None, query_annotation=None,
                 repeat_annotation=None, **kwargs):
        kwargs.update({'use_titles': False})
        super(Ideogram, self).__init__(gff_file=ref_annotation, query=query_annotation,
                                       repeat_annotation=repeat_annotation, **kwargs)
        x_radices, y_radices, x_scale, y_scale = radix_settings  # unpack
        self.border_width = 12
        coordinates = IdeogramCoordinateFrame(x_radices, y_radices, x_scale, y_scale, self.border_width)
        self.all_contents = {}  # (filename: contigs) output_fasta() determines order of fasta_sources
        self.each_layout = [coordinates]  # overwrite anything else
        self.i_layout = 0
        self.layout_algorithm = "1"  # non-raster peano space filling curve
        self.chr_padding = 200


    def process_file(self, input_file_path, output_folder, output_file_name,
                     no_webpage=False, extract_contigs=None):
        self.prep_fasta_source(input_file_path)
        self.image_length = self.calculate_mixed_layout()
        super(Ideogram, self).process_file(input_file_path, output_folder, output_file_name,
                                           no_webpage=no_webpage, extract_contigs=extract_contigs)

    # def activate_high_contrast_colors(self):
    #     # Terrain Colors
    #     self.palette['G'] = hex_to_rgb('6EBAFD')  # Sky or 6EBAFD for darker
    #     self.palette['C'] = hex_to_rgb('EE955D')  # rock
    #     self.palette['T'] = hex_to_rgb('A19E3D')  # light green
    #     self.palette['A'] = hex_to_rgb('6D772F')  # Dark Green

    def render_multiple_files(self, input_fasta_folder, output_folder, output_file_name):
        """Alternative entry point to Ideogram.process_file().  Handles everything."""
        start_time = datetime.now()
        self.final_output_location = output_folder
        make_output_directory(output_folder)
        self.preview_all_files(input_fasta_folder)
        self.image_length = self.calculate_mixed_layout()
        print("Ideogram Point Mapping Done", datetime.now() - start_time)
        self.prepare_image(self.image_length)
        # print("Tallied all contigs :", datetime.now() - start_time)
        print("Initialized Image:", datetime.now() - start_time, "\n")

        for file_no, fasta in enumerate(self.fasta_sources):
            self.i_layout = file_no
            self.contigs = self.all_contents[fasta]  #TODO
            try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
                self.draw_nucleotides()
            except Exception as e:
                print('Encountered exception while drawing nucleotides:', '\n')
                traceback.print_exc()
            try:
                self.draw_extras()
            except BaseException as e:
                print('Encountered exception while drawing extra markup:', '\n')
                traceback.print_exc()

            input_path = os.path.join(input_fasta_folder, fasta)
            self.output_fasta(output_folder, input_path, False, None, False, append_fasta_sources=False)
        print("\nDrew Nucleotides:", datetime.now() - start_time)
        self.output_image(output_folder, output_file_name, False)
        print("Output Image in:", datetime.now() - start_time)


    def draw_nucleotides(self, verbose=True):
        """Only renders one file without modifying the current state of i_layout and self.contigs."""
        total_progess = 0
        for contig in self.contigs:
            if contig.name.lower().startswith('between genes'):
                self.intergenic_palette()
            else:
                self.activate_high_contrast_colors()
            for nuc_i, nucleotide in enumerate(contig.seq):
                try:
                    x, y = self.position_on_screen(total_progess + nuc_i)
                    self.draw_pixel(nucleotide, x, y)
                except IndexError:
                    if verbose:
                        print("Ran out of room at (%i,%i)" % (x,y))
                    break
            total_progess += len(contig.seq)


    def position_on_screen(self, progress):
        return self.levels.position_on_screen(progress)

    def relative_position(self, progress):
        return self.levels.relative_position(progress)

    def draw_extras(self):
        super(Ideogram, self).draw_extras()


    def max_dimensions(self, longest_chr):
        dim = int(np.sqrt(longest_chr * 2))  # ideogram has low density and mostly square
        nucleotide_width = reduce(int.__mul__, self.levels.x_radices)
        y_body = reduce(int.__mul__, self.levels.y_radices[:-1])
        n_coils = np.ceil(longest_chr / nucleotide_width)
        y_needed = int(np.ceil(n_coils / y_body))
        self.levels.y_radices[-1] = y_needed
        n_chr = len(self.each_layout)
        inter_pad = self.chr_padding * (n_chr - 1) if n_chr > 1 else 0
        width = n_chr * (nucleotide_width * self.levels.x_scale + self.levels.origin[0] * 2) + inter_pad
        padding_per_coil = 6
        nuc_height = reduce(int.__mul__, self.levels.y_radices) + padding_per_coil * y_needed
        height = nuc_height * self.levels.y_scale + self.levels.origin[1]*2 + 10

        if self.levels.y_radices[-1] % 2 == 0:  # needs to be odd, but doesn't affect the height
            self.levels.y_radices[-1] += 1
        if self.museum_mode:  # Manually set size for two or more rows of chromosomes
            width, height = 10095, 9820 * 2 + (405*5)
        return width, height


    def draw_extras_for_chromosome(self, scaff_name, coordinate_frame):
        super(Ideogram, self).draw_extras_for_chromosome(scaff_name, coordinate_frame)

    def draw_annotation_labels(self, markup_image, annotated_regions, start_offset, label_color,
                               universal_prefix='', use_suppression=False, force_orientation=None):
        super(Ideogram, self).draw_annotation_labels(markup_image, annotated_regions, start_offset, label_color,
               universal_prefix=universal_prefix, use_suppression=use_suppression,
               force_orientation='horizontal')  # horizontal is important for upper_left

    def draw_label(self, contig_name, width, height, font, title_width, upper_left, vertical_label, strand,
                   canvas, label_color, horizontal_centering=False, center_vertical=False, chop_text=True):
        """Intercept calls from parent and inject some default parameters for Ideograms."""
        # if height < 35 and font.getsize(contig_name)[1] <= 11:
        #     upper_left[1] += 15
        #     upper_left[0] += 15
        self.levels.write_label(contig_name, width, height, font, title_width, upper_left,
              False, '+', canvas, label_color=label_color, horizontal_centering=True, center_vertical=True,
              chop_text=False)

    def calc_all_padding(self):
        """Use the contigs read in to create CoordinateFrames, one for each chromosome."""
        ignored = super(Ideogram, self).calc_all_padding()  # initializes reset values
        # return only the longest chromosome so that it can be used to calculate height
        longest_chr = sum([len(c.seq) for c in self.contigs]) # length used to calculate image size
        return longest_chr

    def calc_padding(self, total_progress, next_segment_length):
        """Ideogram doesn't use within contig padding"""
        return 0,0,0

    def preview_all_files(self, input_fasta_folder):
        """Populates fasta_sources and all_contents with files from a directory"""
        for fasta_file in fastas_in_folder(input_fasta_folder):
            self.prep_fasta_source(fasta_file)
        print("Read all files")

    def prep_fasta_source(self, fasta_file_path):
        self.read_contigs_and_calc_padding(fasta_file_path, None)
        fasta_name = os.path.basename(fasta_file_path)
        self.fasta_sources.append(fasta_name)
        self.all_contents[fasta_name] = self.contigs  # store contigs so they can be wiped

    def calculate_mixed_layout(self):
        from copy import copy
        original_layout = self.each_layout[0]
        origin_series = []
        for chr_count, chromosome in enumerate(self.fasta_sources):
            next_origin = list(original_layout.origin)
            next_origin[0] += (original_layout.base_width + self.chr_padding) * chr_count
            origin_series.append(next_origin)
        if self.museum_mode:  # manual layout by editing code
            #chr1 is 24 coils tall
            #last chr is 7 coils shorter than its neighbor
            coil_height = 405
            origin_series = [(p[0],p[1]) for p in origin_series[:3]]  # add title padding
            origin_series += [(origin_series[0][0], origin_series[0][1]+ coil_height * 27),
                             ((origin_series[1][0], origin_series[1][1]+ coil_height *(27 + 7)))]
        longest_chromosome = 0
        for chr_count, chromosome in enumerate(self.fasta_sources):
            current_length = sum([len(c.seq) for c in self.all_contents[chromosome]])
            longest_chromosome = max(longest_chromosome, current_length)
            coord_frame = copy(original_layout)  # make a shallow copy of coord_frame so we can change origin
            # increment the origin and paste the same coordinateframe as before
            coord_frame.origin = tuple(origin_series[chr_count])
            if chr_count:  # not the first
                self.each_layout.append(coord_frame) # scoot next_origin by width we just used up

        print("Calculating Ideogram Point Mapping")
        # Calculate longest coord_frame first, then reuse
        self.each_layout[0].build_coordinate_mapping(longest_chromosome)
        for layout in self.each_layout:
            layout.point_mapping = self.each_layout[0].point_mapping
        # these should be shallow
        return longest_chromosome

    def intergenic_palette(self):
        self.palette['G'] = hex_to_rgb('FFA889')  # Red
        self.palette['C'] = hex_to_rgb('FFD389')  # Yellow
        self.palette['T'] = hex_to_rgb('8EB1E1')  # Blue originally '0F4FA8'
        self.palette['A'] = hex_to_rgb('89E4B8')  # Green originally ' 00B25C'
        self.palette['N'] = hex_to_rgb('C1C1C1')  # Grey


def increment(digits, radices, place):
    """Manually counting a number where each digit is in a different based determined
    by the corresponding radix number."""
    if digits[place] < (radices[place] - 1):  # still room left in this base
        digits[place] += 1  # increment digit
        return place
    else:
        digits[place] = 0  # hit max value, roll over to next digit on the right
        return increment(digits,radices,place + 1)



if __name__ == "__main__":
    # layout = Ideogram([3,3,3,63], [5,5,3,3,21])
    # layout.process_file("example_data/hg38_chr19_sample.fa", 'www-data/dnadata/test ideogram', 'ideogram-padding2')

    # layout = Ideogram([3,3,3,63], [5,5,3,3,21], 2, 2)
    # layout.process_file("example_data/hg38_chr19_sample.fa", 'www-data/dnadata/test ideogram', 'ideogram-sparse')
    # layout = Ideogram(([5,5,5,5,11],  # thick, local
    #                    [5,5,5,5,5 ,53], 1, 1))
    ## thin layout layout = Ideogram([3,3,3,3,3,27], [3,3,3,3,3,3 ,53], 1, 1)
    #                               ([5,5,5,5,5,27], [3,3,3,3,3,3 ,53], 1, 1)
    # 3*3*3*3*3*27*
    # 3*3*3*3*3 = 1,594,323 bp per fiber row

    radix_settings = eval(sys.argv[2])
    if len(radix_settings) in (2,4):
        assert type(radix_settings[0]) == type(radix_settings[1]) == type([]), \
            "Wrong types: Example: '([5,5,5,5,11], [5,5,5,5,5 ,53])'"
        if len(radix_settings) == 4:
            assert type(radix_settings[2]) == type(radix_settings[3]) == type(1), \
        "Wrong types: Example: '([5,5,5,5,11], [5,5,5,5,5 ,53])'"
    else:
        raise IndexError("Radix should be two lists: Example: '([5,5,5,5,11], [5,5,5,5,5 ,53])'")
    layout = Ideogram(radix_settings)
    input = sys.argv[1]  #r"D:\Genomes\Human\hg38_chr1.fa"  #
    layout.process_file(input,
                        'www-data/dnadata/Ideograms/5,5,5,5,11/',
                        os.path.splitext(input.split('_')[-1])[0])

    beep(200)

    """
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 14" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_14_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 13" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_13_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 12" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_12_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 11" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_11_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 10" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_10_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 09" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_09_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 08" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_08_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 07" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_07_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 06" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_06_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 05" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_05_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 04" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_04_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 03" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_03_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 02" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_02_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --ref_annotation="D:\Genomes\Malaria\Pfalciparum.noseq_filtered.gff3" --outname="Plasmodium falciparum 3D7 01" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],1,1)" --quick --contigs Pf3D7_01_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --outname="Plasmodium falciparum 3D7 API" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],2,2)" --quick --contigs Pf3D7_API_v3
    .\FluentDNA.exe --fasta="D:\Genomes\Malaria\PlasmoDB-39_Pfalciparum3D7_Genome.fasta" --outname="Plasmodium falciparum 3D7 Pf_M76611" --radix="([3,3,3,3, 7], [5,3,3,3,3,53],2,2)" --quick --contigs Pf_M76611
    """