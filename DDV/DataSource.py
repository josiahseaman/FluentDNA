from DNASkittleUtils.Contigs import Contig, write_contigs_to_file, read_contigs
import os
from DDVUtils import write_contigs_to_chunks_dir, filter_by_contigs
from DDV.DDVUtils import multi_line_height, copy_to_sources
from Layouts import LayoutFrame


class DataSource:
    def __init__(self, fasta_name: str, sort_contigs: bool, extract_contigs: bool, layouframe: LayoutFrame):
        self.coords = layouframe  # LayoutFrame
        self.contigs = []  # List[Contig]
        self.fasta_name = fasta_name  # str
        self.spacing_memory = ''
        self.layout_algorithm = "0"  # rastered tile layout
        self.protein_palette = False
        self.using_spectrum = False
        self.skip_small_titles = False
        self.sort_contigs = sort_contigs
        self.extract_contigs = extract_contigs

    def __getitem__(self, index):
        return self.coords[index]

    @property
    def base_width(self):
        """Shorthand for the column width value that is used often.  This can change
        based on the current self.i_layout."""
        return self.coords.base_width

    @property
    def origin(self):
        return self.coords.origin

    def relative_position(self, progress):  # Alias for layout: Optimize?
        return self.coords.relative_position(progress)

    def position_on_screen(self, progress):  # Alias for layout: Optimize?
        return self.coords.position_on_screen(progress)

    def output_fasta(self, output_folder, fasta, no_webpage, extract_contigs,
                     append_fasta_sources=True):
        bare_file = os.path.basename(fasta)

        # also make single file
        if not no_webpage:
            write_contigs_to_chunks_dir(output_folder, bare_file, self.contigs)
            fasta_destination = os.path.join(output_folder, 'sources', bare_file)
            if self.extract_contigs or self.sort_contigs:  # customized_fasta
                length_sum = sum([len(c.seq) for c in self.contigs])
                fasta_destination = '%s__%ibp.fa' % (os.path.splitext(fasta_destination)[0], length_sum)
                write_contigs_to_file(fasta_destination, self.contigs)  # shortened fasta
            else:
                copy_to_sources(output_folder, fasta)
            print("Sequence saved in:", fasta_destination)
        self.clear_sequences()

    def clear_sequences(self):
        self.spacing_memory = self.contig_struct()
        self.contigs = []

    def contig_struct(self):
        if not self.contigs and self.spacing_memory:
            return self.spacing_memory  # original value was already cached
        json = []
        xy_seq_start = 0
        for index, contig in enumerate(self.contigs):
            if index > 1000:
                break  # I don't want to use a slice operator on the for loop because that will copy it
            xy_seq_start += contig.reset_padding + contig.title_padding
            xy_seq_end = xy_seq_start + len(contig.seq)
            json.append(
                {"name": contig.name.replace("'", ""), "xy_seq_start": xy_seq_start, "xy_seq_end": xy_seq_end,
                 "title_padding": contig.title_padding, "tail_padding": contig.tail_padding,
                 "xy_title_start": xy_seq_start - contig.title_padding,
                 "nuc_title_start": contig.nuc_title_start, "nuc_seq_start": contig.nuc_seq_start})
            xy_seq_start += len(contig.seq) + contig.tail_padding
        return json


    def read_contigs_and_calc_padding(self, input_file_path, extract_contigs=None):
        self.extract_contigs = extract_contigs
        self.fasta_name = os.path.basename(input_file_path)
        try:
            self.contigs = read_contigs(input_file_path) # TODO:, extract_contigs)
        except UnicodeDecodeError as e:
            print(e)
            print("Important: Non-standard characters detected.  Switching to 256 colormap for bytes")
            self.using_spectrum = True
            self.contigs = [Contig(input_file_path, open(input_file_path, 'rb').read())]
        self.contigs = filter_by_contigs(self.contigs, extract_contigs)
        self.protein_palette = is_protein_sequence(self.contigs[0])

        # if len(self.levels) >= 5 and len(self.contigs[0].seq) > self.levels[4].chunk_size and multipart_file:
        #     self.enable_fat_headers()  # first contig is huge and there's more contigs coming
        if len(self.contigs) > 10000:
            print("Over 10,000 scaffolds detected!  Titles for entries less than 10,000bp will not be drawn.")
            self.skip_small_titles = True
            # Important! Skipping isn't valid unless they're sorted
            if not self.sort_contigs:
                self.sort_contigs = True
                print("Scaffolds are being sorted by length.")
                # Best to bring the largest contigs to the  forefront
                self.contigs.sort(key=lambda fragment: -len(fragment.seq))


class PaddedContig(Contig):
    def __init__(self, name, seq):
        super(PaddedContig, self).__init__(name, seq)
        self.reset_padding = 0
        self.title_padding = 0
        self.tail_padding = 0


def is_protein_sequence(contig):
    """Checks if there are any peptide characters in the first 100 of the first contig"""
    peptides = {'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'V', 'W', 'X', 'Y'}
    matches = set(contig.seq[:100]).intersection(peptides)
    print("Found matches", matches)
    return len(matches) > 0




