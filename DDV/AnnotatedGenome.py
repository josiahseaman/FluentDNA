import os
from DNASkittleUtils.Contigs import read_contigs

from DDV.Annotations import create_fasta_from_annotation
from DDV.ParallelGenomeLayout import ParallelLayout

class AnnotatedGenomeLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, *args, **kwargs):
        super(AnnotatedGenomeLayout, self).__init__(n_genomes=2, *args, **kwargs)
        self.fasta_file = fasta_file
        self.gff_file = gff_file

    def render_genome(self, output_folder, output_file_name):
        # first just output one scaffold fasta annotaiton
        annotation_fasta = os.path.basename(self.gff_file) + '.fa'
        self.contigs = read_contigs(self.fasta_file)
        chromosomes = [x.name.split()[0] for x in self.contigs]
        lengths = [len(x.seq) for x in self.contigs]
        create_fasta_from_annotation(self.gff_file, chromosomes,
                                     scaffold_lengths=lengths,
                                     out_name=annotation_fasta)

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
