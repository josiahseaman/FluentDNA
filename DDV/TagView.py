from  os.path import join, basename

import math
import numpy as np
from DNASkittleUtils.Contigs import read_contigs

from DDV.Annotations import GFF
from DDV.AnnotatedGenome import AnnotatedGenomeLayout

class TagView(AnnotatedGenomeLayout):
    def __init__(self, fasta_file, gff_file, *args, **kwargs):
        super(TagView, self).__init__(fasta_file, gff_file, *args, **kwargs)
        self.scan_width = 10

    def scan_reverse_complements(self):
        pass

    def output_tag_byte_sequence(self, output_folder, output_file_name):
        assert len(self.contigs) == 1, "Read one sequence first"
        bytes_file = join(output_folder, output_file_name + '__100.bytes')
        seq = self.contigs[0].seq
        n_lines = (len(seq) // self.scan_width)
        match_bytes = bytearray(n_lines * self.scan_width)
        progress = 0
        # calculate oligomer profiles
        # TODO: real code
        oligomer_profiles = np.random.randint(5, size=(n_lines, 3)) #self.scan_width))

        for y, line in enumerate(oligomer_profiles):
            if y + self.scan_width < len(oligomer_profiles):
                for x in range(1, self.scan_width + 1):
                    score = pythonCorrelate(line, oligomer_profiles[y + x])
                    match_bytes[progress] = max(0, int(score * 255))  # don't care about negative values
                    progress += 1
        with open(bytes_file, 'wb') as out:
            out.write(match_bytes)
        return bytes_file

    def render_genome(self, output_folder, output_file_name):
        # empty annotation

        # sequence
        self.contigs = read_contigs(self.fasta_file)

        # tags is viridis
        bytes_file = self.output_tag_byte_sequence(output_folder, output_file_name)

        ### AnnotatedGenome temporary code  ###
        # annotation_fasta = join(output_folder, basename(self.gff_filename) + '.fa')
        # chromosomes = [x.name.split()[0] for x in self.contigs]
        # lengths = [len(x.seq) for x in self.contigs]
        # create_fasta_from_annotation(self.annotation, chromosomes,
        #                              scaffold_lengths=lengths,
        #                              output_path=annotation_fasta)
        super(TagView, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[#annotation_fasta,
                                       self.fasta_file,
                                       bytes_file])

        pass


def pythonCorrelate(x, y):
    n = len(x)
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0.0
    xdiff2 = 0.0
    ydiff2 = 0.0
    for idx in range(n):
        xdiff = float(x[idx]) - avg_x
        ydiff = float(y[idx]) - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    if xdiff2 * ydiff2 == 0.0:
        return 0.0

    return diffprod / math.sqrt(xdiff2 * ydiff2)



def average(values, start=0, length=-1):
    if length == -1:
        length = len(values)
        assert length > 0
        return float(sum(values)) / length
    assert length > 0
    totalSum = 0
    for index in range(start, start + length):
        totalSum += values[index]
    return float(totalSum) / length

