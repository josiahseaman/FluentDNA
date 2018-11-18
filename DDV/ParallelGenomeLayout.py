from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import traceback
from datetime import datetime
from math import floor

from PIL import ImageFont

from DNASkittleUtils.CommandLineUtils import just_the_name
from DDV.TileLayout import TileLayout, font_filename, level_layout_factory


class ParallelLayout(TileLayout):
    def __init__(self, n_genomes, low_contrast=False, base_width=None, column_widths=None):
        # This layout is best used on one chromosome at a time.
        super(ParallelLayout, self).__init__(use_fat_headers=False, sort_contigs=False,
                                             low_contrast=low_contrast, base_width=base_width)
        if column_widths is not None:
            assert len(column_widths) == n_genomes, \
                "Provide the same number of display widths as data sources."
        if column_widths is None:  # just copies the TileLayout levels several times
            column_widths = [self.base_width] * n_genomes

        self.each_layout = []  # one layout per data source assumed same order as self.fasta_sources
        all_columns_height = self.base_width * 10
        columns = self.levels[2]
        cluster_width = sum(column_widths) + columns.padding * n_genomes  # total thickness of data and padding
        cluster_width += columns.padding * 2  # double up on padding between super columns
        column_clusters_per_mega_row = floor(10600 / cluster_width)

        for nth_genome in range(n_genomes):
            standard_modulos = [column_widths[nth_genome], all_columns_height, column_clusters_per_mega_row, 10, 3, 4, 999]
            standard_step_pad = cluster_width - standard_modulos[0] + 6
            standard_padding = [0, 0, standard_step_pad, 6*3, 6*(3**2), 6*(3**3), 6*(3**4)]
            self.each_layout.append(level_layout_factory(standard_modulos, standard_padding))

        self.levels = self.each_layout[0]

        # steps inside a column bundle, not exactly the same as bundles steps
        thicknesses = [self.each_layout[i][0].modulo + 6 for i in range(n_genomes)]
        self.column_offsets = [sum(thicknesses[:i]) for i in range(n_genomes)]

        self.n_genomes = n_genomes
        self.genome_processed = 0
        self.using_background_colors = False
        self.origin = [6, self.levels[3].thickness + 6]  # start with one row for a title, but not subsequent rows
        self.column_colors = "#FFFFFF #E5F3FF #EAFFE5 #FFE7E5 #F8E5FF #FFF3E5 #FFFFE5 #FFF6E5".split()
        self.column_colors = self.column_colors[:self.n_genomes]

    def enable_fat_headers(self):
        pass  # just don't

    def process_file(self, output_folder, output_file_name, fasta_files=list(),
                     no_webpage=False, extract_contigs=None):
        assert len(fasta_files) == self.n_genomes, "List of Genome files must be same length as n_genomes"
        start_time = datetime.now()
        self.image_length = self.read_contigs_and_calc_padding(fasta_files[0], extract_contigs)
        self.prepare_image(self.image_length)
        if self.using_background_colors:
            self.fill_in_colored_borders()
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
        self.draw_the_viz_title(fasta_files)
        self.generate_html(output_folder, output_file_name)  # only furthest right file is downloadable
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)

    def changes_per_genome(self):
        self.levels = self.each_layout[self.genome_processed]
        if self.using_background_colors:
            self.change_background_color(self.genome_processed)

    def position_on_screen(self, progress):
        """ In ParallelLayout, each genome is given a constant x offset in order to interleave the results of each
        genome as it is processed separately.
        """
        x, y = super(ParallelLayout, self).position_on_screen(progress)
        return [x + self.column_offsets[self.genome_processed], y]

    def fill_in_colored_borders(self):
        """When looking at more than one genome, it can get visually confusing as to which column you are looking at.
        To help keep track of it correctly, ParallelGenomeLayout introduces colored borders for each of the columns.
        Then instead of thinking 'I'm looking at the third column' you can think 'I'm looking at the pink column'."""
        # Step through the upper left corner of each column in the file
        column_size = self.levels[2].chunk_size
        margin = 6 // 2
        for genome_index in range(1, self.n_genomes):  # skip the white column
            self.genome_processed = genome_index
            color = self.column_colors[genome_index]
            for column_progress in range(0, self.image_length, column_size):
                left, top = self.position_on_screen(column_progress)
                left, top = max(0, left - margin), max(0, top - margin)
                right, bottom = self.position_on_screen(column_progress + column_size - 1)
                right, bottom = min(self.image.width, right + margin), min(self.image.height, bottom + margin)
                self.draw.rectangle([left, top, right, bottom], fill=color)
        self.genome_processed = 0


    def draw_the_viz_title(self, filenames):
        """Write the names of each of the source files in order so their columns can be identified with their
        column colors"""
        font = ImageFont.truetype(font_filename, 380)
        titles = [just_the_name(x) for x in filenames]  # remove extension and path
        span = '      '.join(titles)
        title_spanning_width = font.getsize(span)[0]  # For centered text
        left_start = self.image.width / 2.0 - title_spanning_width / 2.0
        for genome_index in range(self.n_genomes):
            color = self.column_colors[genome_index]
            title = titles[genome_index]
            text_size = font.getsize(title)
            right = left_start + text_size[0]
            bottom = 6 + text_size[1] * 1.1
            if self.using_background_colors:
                self.draw.rectangle([left_start, 6, right, bottom], fill=color)
            self.draw.text((left_start, 6, right, bottom), title, font=font, fill=(30, 30, 30, 255))
            left_start += font.getsize(title + '      ')[0]


    def change_background_color(self, genome_processed):
        from DDV import gap_char
        def hex_to_rgb(h):
            h = h.lstrip('#')
            return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))

        background = hex_to_rgb(self.column_colors[genome_processed])
        self.palette[gap_char] = background


    def calc_padding(self, total_progress, next_segment_length):
        """Parallel Layouts have a special title which describes the first (main) alignment.
        So padding for their title does not need to be included."""
        # Get original values and level
        reset_padding, title_padding, tail = super(ParallelLayout, self).calc_padding(total_progress,
                                                                                      next_segment_length)
        # no larger than 1 full column or text will overlap
        if title_padding >= self.tile_label_size:
            title_padding = self.levels[2].chunk_size
        # Remove first title
        if total_progress == 0:
            tail += title_padding
            title_padding = 0
            # i = min([i for i in range(len(self.levels)) if next_segment_length + 2600 < self.levels[i].chunk_size])
            # total_padding = total_progress + title_padding + reset_padding + next_segment_length
            # tail = self.levels[i - 1].chunk_size - total_padding % self.levels[i - 1].chunk_size - 1

        return reset_padding, title_padding, tail


    def draw_title(self, total_progress, contig):
        # if total_progress != 0:
        super(ParallelLayout, self).draw_title(total_progress, contig)

    def levels_json(self, ignored):
        """Include only the last layout, with correct origin"""
        return super(ParallelLayout, self).levels_json(self.each_layout[-1])
