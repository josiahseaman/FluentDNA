import traceback
from datetime import datetime

from array import array

from DDVUtils import LayoutLevel, rev_comp, Contig
from RepeatAnnotations import read_repeatmasker_csv, max_consensus_width, blank_line_array, filter_repeats_by_chromosome
from TileLayout import TileLayout


class TransposonLayout(TileLayout):
    def __init__(self):
        super().__init__()
        self.repeat_entries = None


    def create_image_from_preprocessed_alignment(self, input_file_path, consensus_width, num_lines, output_folder, output_file_name):
        self.initialize_image_by_sequence_dimensions(consensus_width, num_lines)  # sets self.layout
        self.read_contigs(input_file_path)
        super(TransposonLayout, self).draw_nucleotides()  # uses self.contigs and self.layout to draw
        self.output_image(output_folder, output_file_name)


    def process_file(self, ref_fasta, output_folder, output_file_name, repeat_annotation_filename=None):
        start_time = datetime.now()
        self.read_all_files(ref_fasta, repeat_annotation_filename)
        consensus_width = max_consensus_width(self.repeat_entries)
        num_lines = len(self.repeat_entries)

        print("Read contigs :", datetime.now() - start_time)

        self.initialize_image_by_sequence_dimensions(consensus_width, num_lines)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)


    def initialize_image_by_sequence_dimensions(self, consensus_width, num_lines):
        self.layout_based_on_repeat_size(consensus_width)
        self.image_length = consensus_width * num_lines
        self.prepare_image(self.image_length)
        print("Image is ", self.image.width, "x", self.image.height)


    def read_all_files(self, ref_fasta, repeat_annotation_filename, column='repName', rep_name='L1PA3'):
        if repeat_annotation_filename is None:  # necessary for inheritance requirements
            raise NotImplementedError("TransposonLayout requires a repeat annotation to work")
        self.repeat_entries = read_repeatmasker_csv(repeat_annotation_filename, column, rep_name)
        self.repeat_entries.sort(key=lambda x: -len(x))  # longest first
        print("Found %i entries under %s" % (len(self.repeat_entries), str(rep_name)))
        self.read_contigs(ref_fasta)


    def layout_based_on_repeat_size(self, consensus_width):
        """change layout to match dimensions of the repeat"""
        self.levels = [
            LayoutLevel("X_in_consensus", consensus_width, 1, 0),  # [0]
            LayoutLevel("Instance_line", consensus_width * 10, consensus_width, 0)  # [1]
        ]
        self.levels.append(LayoutLevel("TypeColumn", 100, padding=20, levels=self.levels))  # [2]
        self.levels.append(LayoutLevel("RowInTile", 10, levels=self.levels))  # [3]
        self.levels.append(LayoutLevel("TileColumn", 3, levels=self.levels))  # [4]
        self.levels.append(LayoutLevel("TileRow", 4, levels=self.levels))  # [5]


    def draw_nucleotides(self):
        processed_contigs = self.create_repeat_fasta_contigs()
        self.contigs = processed_contigs  # TODO: overwriting self.contigs isn't really great data management
        super(TransposonLayout, self).draw_nucleotides()  # uses self.contigs and self.layout to draw


    def create_repeat_fasta_contigs(self):
        consensus_width = max_consensus_width(self.repeat_entries)
        # matches = [contig for contig in self.contigs if contig.name == self.current_chromosome]
        processed_contigs = []
        for contig in self.contigs:
            display_lines = []
            reps_on_chr = filter_repeats_by_chromosome(self.repeat_entries, contig.name)
            for fragment in reps_on_chr:
                line = grab_aligned_repeat(consensus_width, contig, fragment)
                display_lines.append(''.join(line))

            processed_seq = ''.join(display_lines)
            processed_contigs.append(Contig(contig.name, processed_seq, 0, 0, 0,
                                            0, 0))  # TODO: title_length currently doesn't have a title and might break mouse tracking
        return processed_contigs


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