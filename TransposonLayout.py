import traceback
from datetime import datetime
from array import array

from DDVUtils import LayoutLevel, rev_comp, Contig
from RepeatAnnotations import read_repeatmasker_csv, max_consensus_width, blank_line_array
from TileLayout import TileLayout


class TransposonLayout(TileLayout):
    def __init__(self):
        super().__init__()
        self.repeat_entries = None


    def process_file(self, ref_fasta, output_folder, output_file_name, repeat_annotation_filename=None):
        start_time = datetime.now()
        if repeat_annotation_filename is None:  # necessary for inheritance requirements
            raise NotImplementedError("TransposonLayout requires a repeat annotation to work")
        column, rep_name = 'repName', 'L1PA3'
        self.repeat_entries = read_repeatmasker_csv(repeat_annotation_filename, column, rep_name)
        self.repeat_entries.sort(key=lambda x: -len(x))  # longest first
        print("Found %i entries under %s" % (len(self.repeat_entries), str(rep_name)))

        self.layout_based_on_repeat_size()

        self.read_contigs(ref_fasta)
        print("Read contigs :", datetime.now() - start_time)
        self.image_length = max_consensus_width(self.repeat_entries) * len(self.repeat_entries)
        self.prepare_image(self.image_length)
        print("Image is ", self.image.width, "x", self.image.height)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)

    def layout_based_on_repeat_size(self):
        """change layout to match dimensions of the repeat"""
        consensus_width = max_consensus_width(self.repeat_entries)
        self.levels = [
            LayoutLevel("X_in_consensus", consensus_width, 1, 0),  # [0]
            LayoutLevel("Instance_line", consensus_width * 10, consensus_width, 0)  # [1]
        ]
        self.levels.append(LayoutLevel("TypeColumn", 100, padding=20, levels=self.levels))  # [2]
        self.levels.append(LayoutLevel("RowInTile", 10, levels=self.levels))  # [3]
        self.levels.append(LayoutLevel("TileColumn", 3, levels=self.levels))  # [4]
        self.levels.append(LayoutLevel("TileRow", 4, levels=self.levels))  # [5]

    def draw_nucleotides(self):
        consensus_width = max_consensus_width(self.repeat_entries)
        # matches = [contig for contig in self.contigs if contig.name == self.current_chromosome]
        processed_contigs = []

        for contig in self.contigs:
            display_lines = []
            reps_on_chr = [x for x in self.repeat_entries if x.geno_name == contig.name]
            for fragment in reps_on_chr:
                line = blank_line_array(consensus_width, 'X', newline=False)
                nucleotides = fragment.genome_span().sample(contig.seq)
                if fragment.strand == '-':
                    nucleotides = rev_comp(nucleotides)
                if fragment.rep_end - len(nucleotides) < 0:  # sequence I have sampled starts before the beginning of the frame
                    nucleotides = nucleotides[len(nucleotides) - fragment.rep_end:]  # chop off the beginning
                line = line[:fragment.rep_end - len(nucleotides)] + array('u', nucleotides) + line[fragment.rep_end:]
                assert len(line) == consensus_width, "%i, %i, %i" % (len(line), consensus_width, self.repeat_entries.index(fragment))
                display_lines.append(''.join(line))

            processed_seq = ''.join(display_lines)
            processed_contigs.append(Contig(contig.name, processed_seq, 0, 0, 0,
                                            0, 0))  # TODO: title_length currently doesn't have a title and might break mouse tracking

        self.contigs = processed_contigs  # TODO: overwriting self.contigs isn't really great data management
        super(TransposonLayout, self).draw_nucleotides()  # uses self.contigs and self.layout to draw