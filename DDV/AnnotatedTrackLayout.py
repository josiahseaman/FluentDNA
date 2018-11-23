import math
from DNASkittleUtils.Contigs import read_contigs
from itertools import chain
from os.path import join, basename

from DDV.Annotations import create_fasta_from_annotation, GFF, find_universal_prefix, extract_gene_name
from DDV.ParallelGenomeLayout import ParallelLayout
from DDV.HighlightedAnnotation import HighlightedAnnotation
from DDV.DDVUtils import filter_by_contigs


class AnnotatedTrackLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, annotation_width, **kwargs):
        self.annotation_phase = 0  # Means annotations are first, on the left
        columns = [annotation_width, 100]  # TODO: or base_width
        super(AnnotatedTrackLayout, self).__init__(n_genomes=2, column_widths=columns, **kwargs)
        self.fasta_file = fasta_file
        self.gff_filename = gff_file
        self.annotation = GFF(self.gff_filename)

    def render_genome(self, output_folder, output_file_name, extract_contigs=None):
        self.annotation_fasta = join(output_folder, basename(self.gff_filename) +
                                     ('.fa' if extract_contigs is None else '_extracted.fa'))
        self.contigs = read_contigs(self.fasta_file)
        # TODO: Genome is read_contigs twice unnecessarily. This could be sped up.
        self.contigs = filter_by_contigs(self.contigs, extract_contigs)
        extract_contigs = [x.name.split()[0] for x in self.contigs]
        lengths = [len(x.seq) for x in self.contigs]
        create_fasta_from_annotation(self.annotation, extract_contigs,
                                     scaffold_lengths=lengths,
                                     output_path=self.annotation_fasta,
                                     annotation_width=self.annotation_width,
                                     base_width=self.each_layout[1].base_width)
        #check contig filtering
        super(AnnotatedTrackLayout, self).process_file(output_folder,
               output_file_name=output_file_name,
               fasta_files=[self.annotation_fasta, self.fasta_file],
               no_webpage=False, extract_contigs=extract_contigs)

    def changes_per_genome(self):
        super(AnnotatedTrackLayout, self).changes_per_genome()
        self.activate_high_contrast_colors()
        if self.genome_processed == self.annotation_phase:  # Use softer colors for annotations
            self.activate_natural_colors()


    def calc_padding(self, total_progress, next_segment_length):
        """ Skip the exceptions used in Parallel Layouts for first scaffold."""
        reset_padding, title_padding, tail = super(ParallelLayout, self) \
            .calc_padding(total_progress, next_segment_length)
        # no larger than 1 full column or text will overlap
        if title_padding >= self.levels[3].chunk_size:
            title_padding = self.levels[2].chunk_size
        return reset_padding, title_padding, tail


    def draw_extras(self):
        """Drawing Annotations labels"""
        self.i_layout = self.annotation_phase  # restore annotation layout and print labels
        self.genome_processed = self.annotation_phase
        self.read_contigs_and_calc_padding(self.annotation_fasta)
        self.prepare_annotation_labels()



    def draw_the_viz_title(self, fasta_files):
        super(AnnotatedTrackLayout, self).draw_the_viz_title(fasta_files)
        # only draw on the annotation pass, not sequence


    def prepare_annotation_labels(self):
        self.i_layout = self.annotation_phase
        labels = self.annotation.annotations  # dict
        layout = self.contig_struct()
        flattened_annotation = list(chain(*[list(annotation_list) for annotation_list in labels.values()]))
        universal_prefix = find_universal_prefix(flattened_annotation)
        print("Removing Universal Prefix from annotations:", universal_prefix)
        for sc_index, scaffold in enumerate(layout):  # Exact match required (case sensitive)
            scaff_name = scaffold["name"].split()[0]
            if scaff_name not in labels.keys():
                continue
            for entry in labels[scaff_name]:
                if entry.feature in ['gene', 'mRNA']:
                    progress = (entry.start ) // self.base_width *\
                               self.annotation_width + scaffold["xy_seq_start"]
                    end = (entry.end) // self.base_width *\
                               self.annotation_width + scaffold["xy_seq_start"]
                    name = extract_gene_name(entry, universal_prefix)
                    width, height, left, right, top, bottom = \
                        handle_multi_column_annotations(self.levels, progress, end)
                    font_size = 9
                    title_width = 18
                    title_lines = math.ceil(len(name) / title_width)
                    # most gene names aren't unique at 9 characters
                    min_height = 13 if title_lines == 1 else 26
                    height = max(min_height, height)
                    vertical = False #height > width
                    if vertical:
                        width, height = height, width
                        # left, top, right, bottom = bottom, left, top, right
                    font = self.get_font(self.font_name, font_size)
                    self.levels.write_label(name, width, height, font, title_width,
                                            [left, top], vertical, entry.strand, self.image)
        print("Done Drawing annotation labels")

    def additional_html_content(self, html_content):
        """{'CDS':FeatureRep('G', 1),  # 1 priority is the most important
            'exon':FeatureRep('T', 2),
            'gene':FeatureRep('C', 3),
            'mRNA':FeatureRep('A', 4),
            'transcript':FeatureRep('N', 5)}"""
        return {'legend': html_content['legend'] +
                """<p><span><strong>Annotation Colors:</strong>
                <p>Gene = Yellow, mRNA = Green, Exon = Blue, CDS = Red. Gene components are stacked in a hierarchy: 
                 CDS in exons, exons in genes. Only the most exclusive category (CDS) is visible. 
                 Visible yellow regions are introns.  Visible blue (exon, but not CDS) are 3' and 5' UTR.</p></span></p>
                """}  # override in children
    @property
    def annotation_width(self):
        return self.each_layout[self.annotation_phase].levels[0].modulo


def handle_multi_column_annotations(coord_frame, start, stop):
    interval = abs(stop - start)
    upper_left = coord_frame.position_on_screen(start + 2)
    bottom_right = coord_frame.position_on_screen(stop - 2)
    multi_column = abs(bottom_right[0] - upper_left[0]) > coord_frame.base_width
    if True:#multi_column:  # pick the biggest column to contain the label, ignore others
        median_point = interval // 2 + min(start, stop)
        s = median_point // coord_frame.base_width * coord_frame.base_width  # beginning of the line holding median
        left = coord_frame.position_on_screen(s)[0]  # x coordinate of beginning of line
        right = coord_frame.position_on_screen(s + coord_frame.base_width - 2)[0]  # end of one line
        #top is the start or  top of the column
        #bottom is the stop or bottom of the column
        column_step = coord_frame.levels[2].chunk_size
        top_of_column = median_point // column_step * column_step
        top = max(start, top_of_column)
        bottom = min(stop, top_of_column + column_step - 2)
        top, bottom = coord_frame.position_on_screen(top)[1], coord_frame.position_on_screen(bottom)[1]
        height = abs(bottom - top)
    # else:
    #     height = interval // coord_frame.base_width
    width = coord_frame.base_width
    return width, height, left, right, top, bottom
