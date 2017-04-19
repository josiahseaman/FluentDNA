import traceback
from datetime import datetime

from DDVUtils import LayoutLevel
from TileLayout import TileLayout


class TransposonLayout(TileLayout):
    def __init__(self):
        super().__init__()
        self.levels = [
            LayoutLevel("X_in_consensus", 10000, 1, 0),  # [0]
            LayoutLevel("Instance_line", 99999, 10000, 0)  # [1]
        ]

    def process_file(self, input_file_path, output_folder, output_file_name, repeat_annotation_filename=None):
        if repeat_annotation_filename is None:
            raise NotImplementedError("TransposonLayout requires a repeat annotation to work")
        start_time = datetime.now()
        self.image_length = self.read_contigs(input_file_path)
        print("Read contigs :", datetime.now() - start_time)
        self.prepare_image(self.image_length)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        # try:  #not doing more than one repeat type right away
        #     if len(self.contigs) > 1:
        #         print("Drawing %i titles" % len(self.contigs))
        #         self.draw_titles()
        #         print("Drew Titles:", datetime.now() - start_time)
        # except BaseException as e:
        #     print('Encountered exception while drawing titles:', '\n')
        #     traceback.print_exc()
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)


    def draw_nucleotides(self):
        current_line = 0
        repeat_name = "L2"
        for line in self.repeat_entries:  # sorted by chromosome position
            local_point = [(line['repStart']), current_line]
            sequence = self.contigs[0][line['genoStart']: line['genoEnd']]
            self.draw_sequence_line(local_point, self.origin, sequence)