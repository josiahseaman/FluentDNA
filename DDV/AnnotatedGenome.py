from DNASkittleUtils.Contigs import read_contigs

from DDV.Annotations import create_fasta_from_annotation
from DDV.ParallelGenomeLayout import ParallelLayout

class AnnotatedGenomeLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, *args, **kwargs):
        super(AnnotatedGenomeLayout, self).__init__(n_genomes=2, *args, **kwargs)
        self.fasta_file = fasta_file
        self.gff_file = gff_file

    def render_genome(self, output_folder, output_file_name, target_chromosome):
        # first just output one scaffold fasta annotaiton
        out_name = 'Annotation_' + target_chromosome + '.fa'
        create_fasta_from_annotation(self.gff_file,
                                     target_chromosome,
                                     'Annotation_' + target_chromosome + '.fa')

        super(AnnotatedGenomeLayout, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[self.fasta_file, out_name])

    def read_contigs_and_calc_padding(self, input_file_path):
        self.contigs = read_contigs(input_file_path)[:1]

        return self.calc_all_padding()

    def color_changes_per_genome(self):
        if self.genome_processed:  # Use softer colors for annotations
            self.activate_natural_colors()

if __name__ == '__main__':
    layout = AnnotatedGenomeLayout(r"D:\Genomes\Gnetum\Gnetum.final.fa", r"D:\Genomes\Gnetum\Gmm.final.gff")
    layout.render_genome()