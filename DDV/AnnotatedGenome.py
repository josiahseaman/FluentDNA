from  os.path import join, basename

import math
from DNASkittleUtils.Contigs import read_contigs

from DDV.Annotations import create_fasta_from_annotation, GFF
from DDV.ParallelGenomeLayout import ParallelLayout

class AnnotatedGenomeLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, annotation_width=None, *args, **kwargs):
        self.annotation_phase = 0  # Means annotations are first, on the left
        self.annotation_width = annotation_width if annotation_width is not None else 100
        columns = [annotation_width, 100]  # TODO: or base_width
        super(AnnotatedGenomeLayout, self).__init__(n_genomes=2, column_widths=columns, *args, **kwargs)
        self.fasta_file = fasta_file
        self.gff_filename = gff_file
        self.annotation = GFF(self.gff_filename)

    def render_genome(self, output_folder, output_file_name, extract_contigs=None):
        self.annotation_fasta = join(output_folder, basename(self.gff_filename) +
                                     ('.fa' if extract_contigs is None else '_extracted.fa'))
        self.contigs = read_contigs(self.fasta_file)
        # TODO: Genome is read_contigs twice unnecessarily. This could be sped up.
        self.contigs = self.filter_by_contigs(self.contigs, extract_contigs)
        extract_contigs = [x.name.split()[0] for x in self.contigs]
        lengths = [len(x.seq) for x in self.contigs]
        create_fasta_from_annotation(self.annotation, extract_contigs,
                                     scaffold_lengths=lengths,
                                     output_path=self.annotation_fasta,
                                     annotation_width=self.annotation_width,
                                     base_width=self.base_width)
        super(AnnotatedGenomeLayout, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[self.annotation_fasta, self.fasta_file],
                          no_webpage=False, extract_contigs=extract_contigs)


    def read_contigs_and_calc_padding(self, input_file_path, extract_contigs=None):
        self.contigs = read_contigs(input_file_path)
        # TODO: Genome is read_contigs twice unnecessarily. This could be sped up.
        return self.calc_all_padding()


    def changes_per_genome(self):
        self.levels = self.each_layout[self.genome_processed]
        self.activate_high_contrast_colors()
        if self.genome_processed == self.annotation_phase:  # Use softer colors for annotations
            self.activate_natural_colors()


    def calc_padding(self, total_progress, next_segment_length, multipart_file):
        """ Skip the exceptions used in Parallel Layouts for first scaffold."""
        reset_padding, title_padding, tail = super(ParallelLayout, self)\
            .calc_padding(total_progress, next_segment_length, multipart_file)
        # no larger than 1 full column or text will overlap
        if title_padding >= self.levels[3].chunk_size:
            title_padding = self.levels[2].chunk_size
        return reset_padding, title_padding, tail


    def draw_titles(self):
        if self.genome_processed == self.annotation_phase:
            pass
        else:
            super(AnnotatedGenomeLayout, self).draw_titles()  # scaffold names
            self.levels = self.each_layout[self.annotation_phase]
            self.genome_processed = 0
            self.read_contigs_and_calc_padding(self.annotation_fasta)
            self.draw_annotation_labels()



    def draw_the_viz_title(self, fasta_files):
        super(AnnotatedGenomeLayout, self).draw_the_viz_title(fasta_files)
        # only draw on the annotation pass, not sequence
        # self.draw_annotation_labels()

    def draw_annotation_labels(self):
        labels = self.annotation.annotations  # dict
        layout = self.contig_struct()
        for sc_index, scaffold in enumerate(layout):  # Exact match required (case sensitive)
            scaff_name = scaffold["name"].split()[0]
            if scaff_name in labels.keys():
                for entry in labels[scaff_name]:
                    assert isinstance(entry, GFF.Annotation), "This isn't a proper GFF object"
                    if entry.feature == 'gene':
                        progress = (entry.start ) // self.base_width *\
                                   self.annotation_width + scaffold["xy_seq_start"]
                        end = (entry.end) // self.base_width *\
                                   self.annotation_width + scaffold["xy_seq_start"]
                        try:
                            name = entry.attributes['Name']
                        except KeyError:
                            name = ';'.join(['%s=%s' %(key, val) for key, val in entry.attributes.items()])
                        upper_left = self.position_on_screen(progress + 2)
                        bottom_right = self.position_on_screen(end - 2)
                        width = 100  # bottom_right[0] - upper_left[0],
                        font_size = 9
                        title_width = 18
                        title_lines = math.ceil(len(name) / title_width)
                        # most gene names aren't unique at 9 characters
                        min_height = 13 if title_lines == 1 else 26
                        height = max(min_height, bottom_right[1] - upper_left[1])
                        self.write_title(name, width, height, font_size, title_lines, title_width, upper_left,
                                         False, self.image)

        print("Done Drawing annotation labels")

    def additional_html_content(self, html_content):
        return {'legend': html_content['legend'] +
                """<p><span><strong>Annotation Colors:</strong>
                <p>Gene = Blue, Exon = Yellow, CDS = Green. Gene components are stacked in a hierarchy: 
                 CDS in exons, exons in genes. Only the most exclusive category (CDS) is visible. 
                 Visible blue regions are introns.  Visible yellow (exon, but not CDS) are 3' and 5' UTR.</p></span></p>
                """}  # override in children