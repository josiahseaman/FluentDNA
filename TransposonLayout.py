import statistics
import traceback
from collections import defaultdict
from datetime import datetime

import math
from array import array

from DDVUtils import LayoutLevel, rev_comp, Contig
from RepeatAnnotations import read_repeatmasker_csv, max_consensus_width, blank_line_array
from TileLayout import TileLayout


class TransposonLayout(TileLayout):
    def __init__(self):
        super().__init__()
        self.using_mixed_widths = False
        self.repeat_entries = None


    def create_image_from_preprocessed_alignment(self, input_file_path, consensus_width, num_lines, output_folder, output_file_name):
        self.using_mixed_widths = False  # use a consistent consensus_width throughout and standard layout levels
        self.initialize_image_by_sequence_dimensions(consensus_width, num_lines)  # sets self.layout
        self.read_contigs(input_file_path)
        super(TransposonLayout, self).draw_nucleotides()  # uses self.contigs and self.layout to draw
        self.output_image(output_folder, output_file_name)


    def process_file(self, ref_fasta, output_folder, output_file_name, repeat_annotation_filename=None):
        self.using_mixed_widths = True  # we are processing all repeat types with different widths
        start_time = datetime.now()
        self.read_all_files(ref_fasta, repeat_annotation_filename)
        average_width = int(statistics.mean([x.rep_end for x in self.repeat_entries]))  # rough approximation of size
        num_lines = len(self.repeat_entries)
        print("Average Width", average_width, "Entries", num_lines)

        print("Read contigs :", datetime.now() - start_time)

        self.initialize_image_by_sequence_dimensions(average_width, num_lines)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)


    def max_dimensions(self, image_length):
        if self.using_mixed_widths:
            rough = math.ceil(math.sqrt(image_length * 3))
            rough = min(62950, rough)  # hard cap at 4GB images created
            return rough + 50, rough + 50
        else:
            return super(TransposonLayout, self).max_dimensions(image_length)


    def initialize_image_by_sequence_dimensions(self, consensus_width, num_lines):
        self.layout_based_on_repeat_size(consensus_width)
        self.image_length = consensus_width * num_lines
        self.prepare_image(self.image_length)
        print("Image is ", self.image.width, "x", self.image.height)


    def read_all_files(self, ref_fasta, repeat_annotation_filename, column='genoName', rep_name='chr20'):
        if repeat_annotation_filename is None:  # necessary for inheritance requirements
            raise NotImplementedError("TransposonLayout requires a repeat annotation to work")
        self.repeat_entries = read_repeatmasker_csv(repeat_annotation_filename, column, rep_name)
        self.filter_simple_repeats()
        self.repeat_entries.sort(key=lambda x: -len(x) + x.geno_start / 200000000)  # longest first, chromosome position breaks ties
        print("Found %s entries under %s" % ('{:,}'.format(len(self.repeat_entries)), str(rep_name)))
        self.read_contigs(ref_fasta)


    def filter_simple_repeats(self):
        # TODO: option to keep ONLY simple_repeats and delete everything else
        before = len(self.repeat_entries)
        self.repeat_entries = [x for x in self.repeat_entries if x.rep_class != 'Simple_repeat']  # remove 'simple'
        difference = before - len(self.repeat_entries)
        print("Removed", difference, "simple repeats", "{:.1%}".format(difference / before), "of the data.")


    def layout_based_on_repeat_size(self, consensus_width):
        """change layout to match dimensions of the repeat"""
        self.levels = [
            LayoutLevel("X_in_consensus", consensus_width, 1, 0),  # [0]
            LayoutLevel("Instance_line", 400, consensus_width, 0)  # [1]
        ]
        if self.using_mixed_widths:
            self.levels.append(LayoutLevel("TypeColumn", 999999, padding=20, levels=self.levels))  # [2]
            self.levels.append(LayoutLevel("RowInTile", 999999, levels=self.levels))  # [3]
        else:
            self.levels.append(LayoutLevel("TypeColumn", 100, padding=20, levels=self.levels))  # [2]
            self.levels.append(LayoutLevel("RowInTile", 10, levels=self.levels))  # [3]
            self.levels.append(LayoutLevel("TileColumn", 3, levels=self.levels))  # [4]
            self.levels.append(LayoutLevel("TileRow", 4, levels=self.levels))  # [5]


    def draw_nucleotides(self):
        processed_contigs = self.create_repeat_fasta_contigs()
        print("Finished creating contigs")
        self.contigs = processed_contigs  # TODO: overwriting self.contigs isn't really great data management
        self.draw_nucleotides_in_variable_column_width()  # uses self.contigs and self.layout to draw


    def draw_nucleotides_in_variable_column_width(self):
        """Layout a whole set of different repeat types with different widths.  Column height is fixed,
        but column width varies constantly.  Wrapping to the next row is determined by hitting the
        edge of the allocated image."""
        for contig in self.contigs:
            assert contig.consensus_width, "You must set the consensus_width in order to use this layout"
            self.layout_based_on_repeat_size(contig.consensus_width)

            contig_progress = 0
            seq_length = len(contig.seq)
            line_width = contig.consensus_width
            for cx in range(0, seq_length, line_width):
                x, y = self.position_on_screen(contig_progress)
                if x + contig.consensus_width + 10 >= self.image.width:
                    contig_progress = 0  # reset to beginning of line
                    self.skip_to_next_mega_row()
                    x, y = self.position_on_screen(contig_progress)
                if y + self.levels[1].modulo >= self.image.height:
                    print("Ran into bottom of image at", contig.name)
                    return contig  # can't fit anything more
                remaining = min(line_width, seq_length - cx)
                contig_progress += remaining
                for i in range(remaining):
                    nuc = contig.seq[cx + i]
                    # if nuc != 'X':
                    self.draw_pixel(nuc, x + i, y)
            print("Drew", contig.name, "at", self.position_on_screen(contig_progress))
            columns_consumed = math.ceil(contig_progress / self.levels[2].chunk_size)
            self.origin[0] += columns_consumed * self.levels[2].thickness
        print('')


    def skip_to_next_mega_row(self):
        print("Skipping to next row:", self.origin)
        self.origin[0] = self.levels[2].padding  # start at left again
        self.origin[1] += self.levels[3].thickness  # go to next mega row


    def create_repeat_fasta_contigs(self):
        processed_contigs = []
        rep_names = list({x.rep_name for x in self.repeat_entries})
        rep_names.sort()  # iterate through unique set in alphabetical order
        for rep_name in rep_names:
            contig = self.make_contig_from_repName(rep_name)
            lines_in_contig = len(contig.seq) // contig.consensus_width
            # minimum number of repeats based on aspect ratio 1:20
            if lines_in_contig > 10 and lines_in_contig > contig.consensus_width // 20:
                print("Collected repeats sequences for", contig)
                processed_contigs.append(contig)
        return processed_contigs


    def make_contig_from_repName(self, rep_name):
        annotations = [x for x in self.repeat_entries if x.rep_name == rep_name]
        consensus_width = max_consensus_width(annotations)
        ordered_lines = defaultdict(lambda: 'X' * consensus_width)
        for contig in self.contigs:
            for line_number, fragment in enumerate(annotations):
                if fragment.geno_name == contig.name:
                    line = grab_aligned_repeat(consensus_width, contig, fragment)
                    ordered_lines[line_number] = ''.join(line)
        processed_seq = ''.join([ordered_lines[i] for i in range(len(annotations))])
        contig_name = '__'.join([annotations[0].rep_name, annotations[0].rep_family, annotations[0].rep_class])
        return Contig(contig_name, processed_seq, 0, 0, 0,
                      0, 0, consensus_width=consensus_width)


def grab_aligned_repeat(consensus_width, contig, fragment):
    line = blank_line_array(consensus_width, 'X', newline=False)
    nucleotides = fragment.genome_span().sample(contig.seq)
    if fragment.strand == '-':
        nucleotides = rev_comp(nucleotides)
    if fragment.rep_end - len(nucleotides) < 0:  # sequence I have sampled starts before the beginning of the frame
        nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
    line = line[:fragment.rep_end - len(nucleotides)] + array('u', nucleotides) + line[fragment.rep_end:]
    assert len(line) == consensus_width, "%i, %i" % (len(line), consensus_width, )

    return line