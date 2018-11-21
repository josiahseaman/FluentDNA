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

from DNASkittleUtils.Contigs import read_contigs

from DDV.DDVUtils import beep
from DDV.HighlightedAnnotation import HighlightedAnnotation
import os
import numpy as np
from functools import reduce


class Ideogram(HighlightedAnnotation):
    def __init__(self, radix_settings, ref_annotation=None, query_annotation=None,
                 repeat_annotation=None, **kwargs):
        kwargs.update({'use_titles': False})
        super(Ideogram, self).__init__(gff_file=ref_annotation, query=query_annotation,
                                       repeat_annotation=repeat_annotation, **kwargs)
        x_radices, y_radices, x_scale, y_scale = radix_settings  # unpack
        self.x_radices = x_radices
        self.y_radices = y_radices
        self.x_scale, self.y_scale = x_scale, y_scale
        self.point_mapping = [] # for annotation and testing purposes
        self.border_width = 12

    def process_file(self, input_file_path, output_folder, output_file_name,
                     no_webpage=False, extract_contigs=None):
        if extract_contigs is None:
            contigs = read_contigs(input_file_path)
            extract_contigs = [contigs[0].name.split()[0]]
            print("Extracting ", extract_contigs)
        super(Ideogram, self).process_file(input_file_path, output_folder, output_file_name,
                                           no_webpage=no_webpage, extract_contigs=extract_contigs)
    # def activate_high_contrast_colors(self):
    #     # Original DDV Colors
    #     self.palette['G'] = (255, 0, 0)
    #     self.palette['A'] = (0, 255, 0)
    #     self.palette['C'] = (250, 240, 114)
    #     self.palette['T'] = (0, 0, 255)
    #     self.palette['G'] = hex_to_rgb('6EBAFD')  # Sky or 6EBAFD for darker
    #     self.palette['C'] = hex_to_rgb('EE955D')  # rock
    #     self.palette['T'] = hex_to_rgb('A19E3D')  # light green
    #     self.palette['A'] = hex_to_rgb('6D772F')  # Dark Green

    def draw_nucleotides(self):
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
        points_visited = set()
        radices.shape = np.prod(radices.shape)  # flatten

        points_file_name = os.path.join(self.final_output_location, "test_ideogram_points.txt")
        points_file = None # open(points_file_name, 'w')
        if points_file:
            print("Saving locations in {}".format(points_file_name))
        contig = self.contigs[0]  # TODO pluck contig by --contigs
        seq_iter = iter(contig.seq)

        if self.x_scale == 1 and self.y_scale == 1:
            self.draw_loop_optimized(curr_pos, digits, no_pts, parities, radices, seq_iter)
        else:
            prevprev_pos = np.zeros((2,), dtype=np.int) # must start at 0 because of scale multiplaction
            prev_pos = np.zeros((2,), dtype=np.int)  # must start at 0 because of scale multiplaction
            curr_pos = np.zeros((2,), dtype=np.int)
            self.draw_loop_any_scale(curr_pos, digits, no_pts, parities, points_file, points_visited,
                                     prev_pos, prevprev_pos, radices, seq_iter, self.x_scale, self.y_scale)


    def draw_loop_optimized(self, curr_pos, digits, no_pts, parities, radices, seq_iter):
        max_x = reduce(int.__mul__, self.x_radices) - 1 #+ self.origin[0]
        min_x = 0  #self.origin[0]
        odd = 0
        for pts in range(no_pts - 1):
            place = increment(digits, radices, 0)
            parities[0:(place // 2 + 1), place % 2] *= -1
            place += 1
            odd = self.hacked_padding(curr_pos, min_x, max_x, odd, place)
            x = int(curr_pos[1])
            y = int(curr_pos[0])
            curr_pos[place % 2] += parities[place // 2, place % 2]
            try:
                self.draw_pixel(next(seq_iter), x + self.levels.origin[0], y + self.levels.origin[1])
            except StopIteration:
                break  # reached end of sequence
            except IndexError:
                print("Ran out of room at (%i,%i)" % (x,y))
                break
            self.point_mapping.append((x,y))

    def hacked_padding(self, curr_pos, min_x, max_x, odd, place):
        if place % 2 == 0:  # this is an y increments
            if place // 2 == len(self.x_radices) - 1:
                if curr_pos[1] == max_x or curr_pos[1] == min_x:
                    if odd == 1:
                        curr_pos[0] += 3  # y coordinates are in [0]
                    odd = (odd + 1) % 2
        return odd

    def draw_loop_any_scale(self, curr_pos, digits, no_pts, parities, points_file, points_visited, prev_pos,
                            prevprev_pos, radices, seq_iter, x_scale, y_scale):
        for pts in range(no_pts - 1):
            place = increment(digits, radices, 0)
            parities[0:(place // 2 + 1), place % 2] *= -1
            place += 1
            prevprev_pos[:] = prev_pos[:]
            prev_pos[:] = curr_pos[:]
            # assume we move 3 up and 5 across
            x = int(prev_pos[1] * x_scale + self.levels.origin[0])
            y = int(prev_pos[0] * y_scale + self.levels.origin[1])
            if points_file:
                print("{} {}".format(x, y), file=points_file)
            curr_pos[place % 2] += parities[place // 2, place % 2]
            diff = curr_pos - prev_pos
            prev_diff = prev_pos - prevprev_pos
            assert (abs(sum(diff)) == 1)
            # assert (x, y) not in points_visited
            # points_visited.add((x, y))
            self.point_mapping.append((x,y))
            try:
                self.paint_turns(seq_iter, x, y, diff, prev_diff,
                                 prev_pos, prevprev_pos, x_scale, y_scale)
            except IndexError:
                print(x, y, "out of range")
            except StopIteration:
                break  # reached end of sequence


    def paint_turns(self, seq_iter, x, y, diff, prev_diff, prev_pos, prevprev_pos, x_scale, y_scale):
        # right-hand rotation at corner when corner==1, left-hand rotation when corner==1, or no turn (corner == 0)
        turn = prev_diff[0] * diff[1] - prev_diff[1] * diff[0]
        if turn == 0 or (x_scale == 1 and y_scale==1):
            self.draw_pixel(next(seq_iter), x, y)
        if diff[1]:
            # x is changing
            for scale_step in range(1, x_scale):
                self.draw_pixel(next(seq_iter), x + scale_step * int(diff[1]), y)
        elif diff[0]:
            # y is changing
            # NB: underlines will sometimes overwrite previous ones
            x_nudge = int(prevprev_pos[1] - prev_pos[1])
            for scale_step in range(1, y_scale):
                self.draw_pixel(next(seq_iter), x + x_nudge, y + scale_step * int(diff[0]))

    def position_on_screen(self, progress):
        """WARNING: This will not work until after self.draw_loop_optimized
         has populated self.point_mapping"""
        x, y = self.point_mapping[progress]
        return x + self.levels.origin[0], y + self.levels.origin[1]

    def relative_position(self, progress):
        return self.point_mapping[progress]

    def draw_extras(self):
        super(Ideogram, self).draw_extras()


    def max_dimensions(self, image_length):
        dim = int(np.sqrt(image_length * 2))  # ideogram has low density and mostly square
        nucleotide_width = reduce(int.__mul__, self.x_radices)
        y_body = reduce(int.__mul__, self.y_radices[:-1])
        n_coils = np.ceil(image_length / nucleotide_width )
        y_needed = int(np.ceil(n_coils / y_body))
        self.y_radices[-1] = y_needed
        width = nucleotide_width * self.x_scale + self.levels.origin[0] * 2
        padding_per_coil = 6
        nuc_height = reduce(int.__mul__, self.y_radices) + padding_per_coil * y_needed
        height = nuc_height * self.y_scale + self.levels.origin[1]*2 + 10

        if self.y_radices[-1] % 2 == 0:  # needs to be odd, but doesn't affect the height
            self.y_radices[-1] += 1
        return width, height


    def handle_multi_column_annotations(self, region, left, right, top, bottom):
        height = bottom - top
        width = right - left
        return width, height, left, right, top

    def draw_extras_for_chromosome(self, scaff_name, coordinate_frame):
        self.use_titles = True
        super(Ideogram, self).draw_extras_for_chromosome(scaff_name, coordinate_frame)

    def write_label(self, contig_name, width, height, font_size, title_width, upper_left, vertical_label,
                    strand, canvas, horizontal_centering=False, center_vertical=False, chop_text=True,
                    label_color=(50, 50, 50, 255)):
        super(Ideogram, self).write_label(contig_name, width, height, font_size, title_width, upper_left,
                                          False, '+', canvas, horizontal_centering=True, center_vertical=True,
                                          chop_text=False)

    def levels_json(self, ignored):
        return '[]'  # There's no reasonable way to encode mouse position in rectangles
    def contig_json(self):
        return '[]'  # There's no reasonable way to encode mouse position in rectangles


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
    assert len(radix_settings) == 4 and \
            type(radix_settings[0]) == type(radix_settings[1]) == type([]) and \
            type(radix_settings[2]) == type(radix_settings[3]) == type(1), \
        "Wrong types: Example: '([5,5,5,5,11], [5,5,5,5,5 ,53], 1, 1)'"
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