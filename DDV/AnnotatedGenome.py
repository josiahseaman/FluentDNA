from  os.path import join, basename

import math
from DNASkittleUtils.Contigs import read_contigs

from DDV.Annotations import create_fasta_from_annotation, GFF
from DDV.ParallelGenomeLayout import ParallelLayout

class AnnotatedGenomeLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, *args, **kwargs):
        super(AnnotatedGenomeLayout, self).__init__(n_genomes=2, *args, **kwargs)
        self.fasta_file = fasta_file
        self.gff_filename = gff_file
        self.annotation = GFF(self.gff_filename) if self.gff_filename is not None else None

    def render_genome(self, output_folder, output_file_name):
        annotation_fasta = join(output_folder, basename(self.gff_filename) + '.fa')
        self.contigs = read_contigs(self.fasta_file)
        chromosomes = [x.name.split()[0] for x in self.contigs]
        lengths = [len(x.seq) for x in self.contigs]
        create_fasta_from_annotation(self.annotation, chromosomes,
                                     scaffold_lengths=lengths,
                                     output_path=annotation_fasta)
        super(AnnotatedGenomeLayout, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[annotation_fasta, self.fasta_file])

    #
    # def read_contigs_and_calc_padding(self, input_file_path):
    #     self.contigs = read_contigs(input_file_path)
    #     # TODO: Genome is read_contigs twice unnecessarily. This could be sped up.
    #     return self.calc_all_padding()


    def color_changes_per_genome(self):
        if not self.using_spectrum:
            self.activate_high_contrast_colors()
            if not self.genome_processed:  # Use softer colors for annotations
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
        super(AnnotatedGenomeLayout, self).draw_titles()  # scaffold names
        if not self.genome_processed:  # only draw on the annotation fasta pass, not sequence
            self.draw_annotation_labels()

    def draw_annotation_labels(self):
        labels = self.annotation.annotations  # dict
        layout = self.contig_struct()
        for sc_index, scaffold in enumerate(layout):  # Exact match required (case sensitive)
            scaff_name = scaffold["name"]
            if scaff_name in labels.keys():
                for entry in labels[scaff_name]:
                    assert isinstance(entry, GFF.Annotation), "This isn't a proper GFF object"
                    if entry.feature == 'gene':
                        progress = (int(entry.start) // 100) * 100 + 2 + scaffold["xy_seq_start"]
                        end = int(entry.end) + scaffold["xy_seq_start"]
                        try:
                            name = entry.attributes['Name']
                        except KeyError:
                            name = ';'.join(['%s=%s' %(key, val) for key, val in entry.attributes.items()])
                        upper_left = self.position_on_screen(progress)
                        bottom_right = self.position_on_screen(end - 2)
                        width = 100  # bottom_right[0] - upper_left[0],
                        font_size = 9
                        title_width = 18
                        title_lines = math.ceil(len(name) / title_width)
                        # most gene names aren't unique at 9 characters
                        min_height = 13 if title_lines == 1 else 26
                        height = max(min_height, bottom_right[1] - upper_left[1])
                        self.write_title(name, width, height, font_size, title_lines, title_width, upper_left, False)

        print("Done Drawing annotation labels")

