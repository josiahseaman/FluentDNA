from  os.path import join, basename
from DNASkittleUtils.Contigs import read_contigs

from DDV.Annotations import create_fasta_from_annotation
from DDV.ParallelGenomeLayout import ParallelLayout

class AnnotatedGenomeLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, *args, **kwargs):
        super(AnnotatedGenomeLayout, self).__init__(n_genomes=2, *args, **kwargs)
        self.fasta_file = fasta_file
        self.gff_file = gff_file

    def render_genome(self, output_folder, output_file_name):
        annotation_fasta = join(output_folder, basename(self.gff_file) + '.fa')
        self.contigs = read_contigs(self.fasta_file)
        chromosomes = [x.name.split()[0] for x in self.contigs]
        lengths = [len(x.seq) for x in self.contigs]
        create_fasta_from_annotation(self.gff_file, chromosomes,
                                     scaffold_lengths=lengths,
                                     output_path=annotation_fasta)

        super(AnnotatedGenomeLayout, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[annotation_fasta, self.fasta_file])

    def read_contigs_and_calc_padding(self, input_file_path):
        self.contigs = read_contigs(input_file_path)
        # TODO: Genome is read_contigs twice unnecessarily. This could be sped up.
        return self.calc_all_padding()

    def color_changes_per_genome(self):
        self.activate_high_contrast_colors()
        if not self.genome_processed:  # Use softer colors for annotations
            self.activate_natural_colors()

    def calc_padding(self, total_progress, next_segment_length, multipart_file):
        """ Skip the exceptions used in Parallel Layouts for first scaffold."""
        reset_padding, title_padding, tail = super(ParallelLayout, self)\
            .calc_padding(total_progress, next_segment_length, multipart_file)
        # no larger than 1 full column or text will overlap
        if title_padding >= self.tile_label_size:
            title_padding = self.levels[2].chunk_size
        return reset_padding, title_padding, tail


    def draw_labels(self, labels_filename):
        import json
        if labels_filename is None:
            return

        labels = json.loads(open(labels_filename, 'r').read())
        for gene in labels:
            start = (int(gene['start']) // 100) * 100 + 2
            end = int(gene['end'])
            name = gene['gene_name']
            upper_left = self.position_on_screen(start)
            bottom_right = self.position_on_screen(end - 2)
            height = max(10, bottom_right[1] - upper_left[1])
            width = 100  # bottom_right[0] - upper_left[0],
            height = min(100, max(22, bottom_right[1] - upper_left[1]))
            font_size = 18
            title_width = 9
            title_lines = 2  # math.ceil(len(name) / title_width)
            self.write_title(name, width, height, font_size, title_lines, title_width, upper_left, False)

        print("Done Drawing annotation labels")

