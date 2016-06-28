from datetime import datetime
from math import floor

from os import path

from LargeImages import DDVTileLayout, LayoutLevel


class ParallelLayout(DDVTileLayout):
    def __init__(self, n_genomes):
        super(ParallelLayout, self).__init__()
        self.use_fat_headers = False  # This layout is best used on one chromosome at a time.
        self.n_genomes = n_genomes
        self.genome_processed = 0
        # modify layout with an additional bundled column layer
        columns = self.levels[2]
        new_width = columns.thickness * n_genomes + columns.padding * 2
        self.levels = self.levels[:2]  # trim off the others
        self.levels.append(LayoutLevel("ColumnInRow", floor(10600 / new_width), levels=self.levels))  # [2]
        self.levels[2].thickness = new_width  # 100+6+100+6+100+18 = 330  #total row width of 10,560 vs original 10,600
        self.levels[2].padding = new_width - (columns.thickness - columns.padding)
        self.column_offset = columns.thickness  # steps inside a column bundle, not exactly the same as bundles steps
        # because of inter bundle padding of 18 pixels
        self.levels.append(LayoutLevel("RowInTile", 10, levels=self.levels))  # [3]
        self.levels.append(LayoutLevel("XInTile", 3, levels=self.levels))  # [4]
        self.levels.append(LayoutLevel("YInTile", 99, levels=self.levels))  # [5]


    def process_file(self, file1, output_folder, output_file_name, additional_files=[]):
        assert len(additional_files) + 1 == self.n_genomes, "List of Genome files must be same length as n_genomes"
        start_time = datetime.now()
        self.image_length = self.read_contigs(file1)
        self.image_length = max(self.image_length, *[path.getsize(file) for file in additional_files])
        print("Read first sequence :", datetime.now() - start_time)
        self.prepare_image(self.image_length)
        self.fill_in_colored_borders()
        print("Initialized Image:", datetime.now() - start_time)
        self.draw_nucleotides()
        print("Drew First File:", file1, datetime.now() - start_time)

        try:
            # Do inner work for two other files
            for filename in additional_files:
                self.genome_processed += 1
                self.read_contigs(filename)
                self.draw_nucleotides()
                print("Drew Additional File:", filename, datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n', str(e))
        self.generate_html(file1, output_folder, output_file_name)
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)


    def position_on_screen(self, index):
        """ In ParallelLayout, each genome is given a constant x offset in order to interleave the results of each
        genome as it is processed separately.
        """
        x, y = super(ParallelLayout, self).position_on_screen(index)
        return [x + self.column_offset * self.genome_processed, y]

    def draw_titles(self):
        return  # Titles should not be draw in 3 column layout

    def calc_padding(self, total_progress, next_segment_length, multipart_file):
        """Shallow copy of super().calc_padding where space is not allocated for titles"""
        return 0, 0, 0
        # min_gap = (20 + 6) * 100  # 20px font height, + 6px vertical padding  * 100 nt per line
        # if not multipart_file:
        #     return 0, 0, 0
        #
        # for i, current_level in enumerate(self.levels):
        #     if next_segment_length < current_level.chunk_size:
        #         # give a full level of blank space just in case the previous
        #         title_padding = 0 if self.levels[i - 1].chunk_size > min_gap else min_gap
        #         space_remaining = current_level.chunk_size - total_progress % current_level.chunk_size
        #         # sequence comes right up to the edge.  There should always be >= 1 full gap
        #         reset_level = current_level  # bigger reset when close to filling chunk_size
        #         if next_segment_length + title_padding < space_remaining:
        #             reset_level = self.levels[i - 1]
        #             # reset_level = self.levels[i - 1]
        #         reset_padding = 0
        #         if total_progress != 0:  # fill out the remainder so we can start at the beginning
        #             reset_padding = reset_level.chunk_size - total_progress % reset_level.chunk_size
        #         total_padding = total_progress + title_padding + reset_padding + next_segment_length
        #         tail = self.levels[i - 1].chunk_size - total_padding % self.levels[i - 1].chunk_size - 1
        #
        #         return reset_padding, title_padding, tail

        # return 0, 0, 0

    def fill_in_colored_borders(self):
        """When looking at more than one genome, it can get visually confusing as to which column you are looking at.
        To help keep track of it correctly, ParallelGenomeLayout introduces colored borders for each of the columns.
        Then instead of thinking 'I'm looking at the third column' you can think 'I'm looking at the pink column'."""
        column_colors = "#FFF #fbb4ae #b3cde3 #ccebc5 #decbe4 #fed9a6 #ffffcc #e5d8bd".split()
        column_colors = column_colors[:self.n_genomes]
        # Step through the upper left corner of each column in the file
        column_size = self.levels[2].chunk_size
        self.genome_processed = 0
        for genome_index in range(self.n_genomes):
            color = column_colors[genome_index]
            for column_progress in range(0, self.image_length, column_size):
                left, top = self.position_on_screen(column_progress)
                left, top = max(0, left - 3), max(0, top - 3)
                right, bottom = self.position_on_screen(column_progress + column_size)
                right, bottom = min(self.image.width, right + 3 + 1), min(self.image.height, bottom + 3 + 1)
                self.draw.rectangle([left, top, right, bottom], fill=color)
            self.genome_processed += 1
        self.genome_processed = 0