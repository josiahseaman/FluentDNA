import math
from DNASkittleUtils.Contigs import read_contigs
from itertools import chain
from os.path import join, basename

from DDV.Annotations import create_fasta_from_annotation, find_universal_prefix, parseGFF
from DDV.ParallelGenomeLayout import ParallelLayout
from DDV.DDVUtils import filter_by_contigs, copy_to_sources


class AnnotatedTrackLayout(ParallelLayout):
    def __init__(self, fasta_file, gff_file, annotation_width, **kwargs):
        self.annotation_phase = 0  # Means annotations are first, on the left
        self.genome_phase = 1  #Genome is second, on the right
        columns = [annotation_width, 100]  # TODO: or base_width
        super(AnnotatedTrackLayout, self).__init__(n_genomes=2, column_widths=columns, **kwargs)
        self.fasta_file = fasta_file
        self.gff_filename = gff_file
        self.annotation = parseGFF(self.gff_filename)

    def render_genome(self, output_folder, output_file_name, extract_contigs=None):
        self.annotation_fasta = join(output_folder, 'sources', basename(self.gff_filename) +
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
                                     base_width=self.each_layout[self.genome_phase].base_width)
        #check contig filtering
        super(AnnotatedTrackLayout, self).process_file(output_folder,
               output_file_name=output_file_name,
               fasta_files=[self.annotation_fasta, self.fasta_file],
               no_webpage=False, extract_contigs=extract_contigs)
        # save original GFF for reproducibility
        copy_to_sources(output_folder, self.gff_filename)

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
        genome_width = self.each_layout[self.genome_phase].base_width
        self.i_layout = self.annotation_phase
        labels = self.annotation  # dict
        layout = self.contig_struct()
        genes_seen = set()
        flattened_annotation = list(chain(*[list(annotation_list) for annotation_list in labels.values()]))
        universal_prefix = find_universal_prefix(flattened_annotation)
        print("Removing Universal Prefix from annotations:", universal_prefix)
        for sc_index, scaffold in enumerate(layout):  # Exact match required (case sensitive)
            scaff_name = scaffold["name"].split()[0]
            if scaff_name not in labels.keys():
                continue
            for entry in labels[scaff_name]:
                if entry.type in ['gene', 'mRNA'] and \
                        (not entry.parent() or entry.parent() not in genes_seen): #mRNA double of a gene
                    genes_seen.add(entry.id())
                    name = entry.name(universal_prefix)
                    progress = (entry.start ) // genome_width *\
                               self.annotation_width + scaffold["xy_seq_start"]
                    end = (entry.end) // genome_width *\
                               self.annotation_width + scaffold["xy_seq_start"]
                    if name == '989535g01':
                        print(name, progress)
                    width, height, left, right, top, bottom = \
                        self.levels.handle_multi_column_annotations(progress, end)
                    font_size = 9
                    title_width = 18
                    title_lines = math.ceil(len(name) / title_width)
                    # most gene names aren't unique at 9 characters
                    min_height = 13 if title_lines == 1 else 26
                    height = max(min_height, height)
                    vertical = height > width
                    if vertical:
                        width, height = height, width
                        # left, top, right, bottom = bottom, left, top, right
                    old_with = width
                    width = max(title_width * 6, width)  # Don't truncate the width such that no meaningful text shows up
                    if vertical and entry.strand == '-':  # Anchored on uppef_left, affected by rotation
                        top -= max(0, abs(width - old_with))
                    font = self.get_font(font_size)
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

    def legend(self):
        return "<strong>Legend:</strong>" + \
                self.legend_line('Adenine (A)', 'A') +\
                self.legend_line('Thymine (T)', 'T') +\
                self.legend_line('Guanine (G)', 'G') +\
                self.legend_line('Cytosine (C)', 'C') +\
                self.legend_line('Unsequenced', 'N') +\
                "<strong>Annotation Legend:</strong>"+\
                self.legend_line('Aspartic acid (D)', 'D') +\
                self.legend_line('Glutamic acid (E)', 'E') +\
                self.legend_line('Phenylalanine (F)', 'F') +\
                self.legend_line('Histidine (H)', 'H') +\
                self.legend_line('Isoleucine (I)', 'I') +\
                self.legend_line('Lysine (K)', 'K') +\
                self.legend_line('Leucine (L)', 'L') +\
                self.legend_line('Methionine (M)', 'M') +\
                self.legend_line('Asparagine (N)', 'N') +\
                self.legend_line('Proline (P)', 'P') +\
                self.legend_line('Glutamine (Q)', 'Q') +\
                self.legend_line('Arginine (R)', 'R') +\
                self.legend_line('Serine (S)', 'S') +\
                self.legend_line('Valine (V)', 'V') +\
                self.legend_line('Tryptophan (W)', 'W') +\
                self.legend_line('Tyrosine (Y)', 'Y')+ \
                self.legend_line('Any (X)', 'X')
