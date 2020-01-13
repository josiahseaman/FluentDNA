from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import os
import traceback
from datetime import datetime

from PIL import ImageFont, Image

from DNASkittleUtils.CommandLineUtils import just_the_name
from FluentDNA.TileLayout import TileLayout, hex_to_rgb
from FluentDNA.Layouts import level_layout_factory


class ParallelLayout(TileLayout):
    def __init__(self, n_genomes, low_contrast=False, base_width=100, column_widths=None,
                 border_boxes=False):
        # This layout is best used on one chromosome at a time.
        super(ParallelLayout, self).__init__(sort_contigs=False,
                                             low_contrast=low_contrast, base_width=base_width)
        self.use_border_boxes = border_boxes
        self.header_height = 12 if border_boxes else 0

        if column_widths is not None:
            assert len(column_widths) == n_genomes, \
                "Provide the same number of display widths as data sources."
        if column_widths is None:  # just copies the TileLayout levels several times
            column_widths = [self.base_width] * n_genomes

        self.each_layout = []  # one layout per data source assumed same order as self.fasta_sources
        # I found that less padding is better for keeping visual patterns coherent over clusters
        # of columns.  The white space has a disproportionate effect if you space it out too much.
        p = 1 if not self.use_border_boxes else 6  # padding_between_layouts
        cluster_width = sum(column_widths) + p * n_genomes  # total thickness of data and padding
        cluster_width += p * 2  # double up on padding between super columns
        column_clusters_per_mega_row = 10600 // cluster_width  # 10600

        for nth_genome in range(n_genomes):
            all_columns_height = base_width * 10
            standard_modulos = [column_widths[nth_genome], all_columns_height, column_clusters_per_mega_row, 109, 999]
            standard_step_pad = cluster_width - standard_modulos[0] + p
            mega_row_padding = p * 3 + self.header_height + 10
            standard_padding = [0, 0, standard_step_pad, mega_row_padding, 777]
            # steps inside a column bundle, not exactly the same as bundles steps
            thicknesses = [other_layout[0].modulo + p for other_layout in self.each_layout]
            origin = (sum(thicknesses) + p, p + self.header_height)
            self.each_layout.append(level_layout_factory(standard_modulos, standard_padding, origin))

        self.n_genomes = n_genomes
        self.genome_processed = 0
        self.megarow_label_size = self.levels[3].chunk_size

    def process_file(self, output_folder, output_file_name, fasta_files,
                     no_webpage=False, extract_contigs=None):
        assert len(fasta_files) == self.n_genomes, "List of Genome files must be same length as n_genomes"
        start_time = datetime.now()
        self.image_length = self.read_contigs_and_calc_padding(fasta_files[0], extract_contigs)
        self.prepare_image(self.image_length)
        if self.use_border_boxes:
            self.draw_border_boxes(fasta_files)
        print("Initialized Image:", datetime.now() - start_time)

        try:
            # Do inner work for each file
            for index, filename in enumerate(fasta_files):
                self.changes_per_genome()
                if index != 0:
                    self.read_contigs_and_calc_padding(filename, extract_contigs)
                self.draw_nucleotides()
                if index == self.n_genomes -1: #last one
                    self.draw_titles()
                self.genome_processed += 1
                print("Drew File:", filename, datetime.now() - start_time)
                self.output_fasta(output_folder, filename, False, extract_contigs, self.sort_contigs)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        try:
            self.draw_extras()
        except BaseException as e:
            print('Encountered exception while drawing titles:', '\n')
            traceback.print_exc()
        # self.draw_the_viz_title(fasta_files)  # Needs padding in origins to work
        # self.generate_html(output_folder, output_file_name) # done in fluentdna.py
        self.output_image(output_folder, output_file_name, no_webpage)
        print("Output Image in:", datetime.now() - start_time)
        return start_time

    def changes_per_genome(self):
        self.i_layout = self.genome_processed

    def position_on_screen(self, progress):
        """ In ParallelLayout, each genome is given a constant x offset in order to interleave the results of each
        genome as it is processed separately.
        """
        x, y = super(ParallelLayout, self).position_on_screen(progress)
        return [x, y]

    def draw_border_boxes(self, fasta_files):
        """When looking at more than one genome, it can get visually confusing as to which column you are looking at.
        To help keep track of it correctly, ParallelGenomeLayout demarcates bundles of columns that go
        together.  Mouse over gives further information on each file."""
        from DNASkittleUtils.DDVUtils import pp
        from FluentDNA.FluentDNAUtils import execution_dir
        base_dir = execution_dir()
        # Caution: These corners are currently hard coded to the color and dimension of one image
        try:
            corner = Image.open(os.path.join(base_dir,'FluentDNA','html_template','img','border_box_corner.png'))
        except (FileNotFoundError, NotADirectoryError):
            try:
                corner = Image.open(os.path.join(base_dir, 'html_template', 'img', 'border_box_corner.png'))
            except (FileNotFoundError, NotADirectoryError):
                corner = Image.open(os.path.join(os.path.dirname(base_dir),
                                                 'html_template', 'img', 'border_box_corner.png'))
        corner_rb = corner.copy().rotate(270, expand=True)
        corner_lb = corner.copy().rotate(180, expand=True)
        corner_lt = corner.copy().rotate(90, expand=True)
        column_size = self.levels[2].chunk_size
        margin = 6
        color = hex_to_rgb('#c9c9c9')
        main_contig = self.contigs[0]
        for column_progress in range(main_contig.title_padding + main_contig.reset_padding,
                                     len(main_contig.seq) + main_contig.title_padding + column_size, column_size):
            left, top = self.each_layout[0].position_on_screen(column_progress)
            left, top = max(0, left - margin), max(0, top - margin - self.header_height)
            # column_progress only works when first and last columns have the same width
            last_column = self.each_layout[-1]
            right, bottom = last_column.position_on_screen(column_progress + column_size - 1)
            right, bottom = min(self.image.width, right + margin), min(self.image.height, bottom + margin//2)
            self.draw.rectangle([left, top, right, bottom], fill=color)
            self.image.paste(corner, (right -6, top))
            self.image.paste(corner_rb, (right - 7, bottom - 6))
            self.image.paste(corner_lb, (left , bottom - 7))
            self.image.paste(corner_lt, (left, top))
            #TODO: could be optimized by caching the text image
            for i, layout in enumerate(self.each_layout):
                left, ignore = layout.position_on_screen(column_progress)
                text = pp(column_progress - main_contig.title_padding) #+ ' ' + just_the_name(fasta_files[i]) # cluttered
                self.write_title(text, self.base_width, self.header_height + 2, 11, 1, 30,
                                 (left, top + 0),
                                 False, self.image, color=hex_to_rgb('#606060'))
        self.genome_processed = 0


    def draw_the_viz_title(self, filenames):
        """Write the names of each of the source files in order so their columns can be identified with their
        column colors"""
        font = self.get_font(380)
        titles = [just_the_name(x) for x in filenames]  # remove extension and path
        span = '      '.join(titles)
        title_spanning_width = font.getsize(span)[0]  # For centered text
        left_start = self.image.width / 2.0 - title_spanning_width / 2.0
        for genome_index in range(self.n_genomes):
            title = titles[genome_index]
            text_size = font.getsize(title)
            right = left_start + text_size[0]
            bottom = 6 + text_size[1] * 1.1
            self.draw.text((left_start, 6, right, bottom), title, font=font, fill=(30, 30, 30, 255))
            left_start += font.getsize(title + '      ')[0]


    def calc_padding(self, total_progress, next_segment_length):
        """Parallel Layouts have a special title which describes the first (main) alignment.
        So padding for their title does not need to be included."""
        # Get original values and level
        reset_padding, title_padding, tail = super(ParallelLayout, self).calc_padding(total_progress,
                                                                                      next_segment_length)
        # no larger than 1 full column or text will overlap
        if title_padding >= self.levels[2].chunk_size:
            title_padding = self.levels[2].chunk_size
        # Remove first title
        # if total_progress == 0:
        #     # tail += title_padding  #commenting this out could cause problems with multiple contigs
        #     title_padding = 0
        #     if len(self.contigs) == 1:
        #         tail = 0  # no need for tail
        #     i = min([i for i in range(len(self.levels)) if next_segment_length + 2600 < self.levels[i].chunk_size])
        #     total_padding = total_progress + title_padding + reset_padding + next_segment_length
        #     tail = self.levels[i - 1].chunk_size - total_padding % self.levels[i - 1].chunk_size - 1

        return reset_padding, title_padding, tail


    def draw_title(self, total_progress, contig):
        super(ParallelLayout, self).draw_title(total_progress, contig)
