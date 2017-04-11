import argparse
import math
from collections import defaultdict
import sys

from DDVUtils import make_output_dir_with_suffix, base_directories
from DDV import create_tile_layout_viz_from_fasta
from TileLayout import TileLayout


def hasDepth(listLike):
    try:
        return len(listLike) > 0 and not isinstance(listLike, (str, dict, tuple, type(u"unicode string"))) and hasattr(
            listLike[0], "__getitem__")
    except:
        return False


def interpolate(A, B, start, end, position):
    if start == end:
        return A
    progress = (position - start) / (end - start)  # progress goes from 0.0 p1  to 1.0 p2
    inverse = 1.0 - progress
    sample = A * inverse + B * progress
    return sample


def linspace(start, end, steps):
    return [interpolate(start, end, 0, steps - 1, i) for i in range(steps)]


def countNucleotides(seq, oligomerSize):
    if hasDepth(seq):
        return [countNucleotides(x, oligomerSize) for x in seq if x != '' and x != [] and x != {}]
    if not seq:
        return {}
    counts = defaultdict(lambda: 0)
    for endIndex in range(oligomerSize, len(seq) + 1, 1):
        c = seq[endIndex - oligomerSize: endIndex]
        counts[c] += 1  # defaults to 0
    return counts


class Sequenaut(TileLayout):  # TODO: make an abstract class parent
    def __init__(self, oligomer_size=100, peak=10.0, baseline=1.0):
        super(Sequenaut, self).__init__()
        self.oligomer_size = oligomer_size
        self.peak = peak
        self.baseline = baseline
        midpoint = self.oligomer_size // 2
        self.weight_matrix = list(linspace(1.0, self.peak, midpoint)) + list(linspace(self.peak, 1.0, self.oligomer_size - midpoint))
        self.image_length = self.max_dimensions(oligomer_size, verbose=True)[0]
        self.hits = [[0] * self.image_length for i in range(self.image_length)]

    @staticmethod
    def coordinate(olig, weight_matrix):
        """Turns an Olig into an (x,y) tuple coordinate.
        axes = {'A': Point(0,1), 'G': Point(0,0), 'C': Point(1,0), 'T': Point(1,1)}"""
        x = 0.0
        y = 0.0
        for c, weight in zip(olig, weight_matrix):
            if c == 'T':
                x += weight
                y += weight
            elif c == 'A':
                y += weight
            elif c == 'C':
                x += weight
        return int(x), int(y)


    def max_dimensions(self, ignored, verbose=False):
        resolution = self.coordinate('T' * self.oligomer_size, self.weight_matrix)[0] + 1
        if verbose:
            print("Image will be %i x %i: %s pixels" % (resolution, resolution, "{:,}".format(resolution**2)))
        return [resolution, resolution]


    def weighted_sequenaut(self, seq):
        counts = countNucleotides(seq, self.oligomer_size)
        for olig, count in counts.items():
            point = self.coordinate(olig, self.weight_matrix)
            self.hits[point[1]][point[0]] += count


    def draw_nucleotides(self):
        for contig in self.contigs:
            self.weighted_sequenaut(contig.seq)
        leader = math.log(max([max(column) for column in self.hits]))
        try:
            middle = self.max_dimensions(self.oligomer_size)[0] // 2
            print("Leader", leader, int(math.log(self.hits[middle][middle]) / leader * 235 + 20))
        except: pass  # just talking

        for y in range(len(self.hits)):
            for x in range(len(self.hits[y])):
                if self.hits[y][x]:
                    grey = 255 - int(math.log(self.hits[y][x]) / leader * 235 + 20)
                    self.pixels[x, y] = (grey, grey, grey)


    def draw_titles(self):
        pass  # There are no titles in Sequenaut


    def process_file(self, input_file_path, output_folder, output_file_name):
        # use the parent process_file() but override methods in child
        super(Sequenaut, self).process_file(input_file_path, output_folder, output_file_name)


def run_sequenaut(args):
    SERVER_HOME, base_path = base_directories(args)
    if args.layout_type == "triangle_slide":  # Typical Use Case
        # TODO: allow batch of tiling layout by chromosome
        output_dir = make_output_dir_with_suffix(base_path, '')
        renderer = Sequenaut(oligomer_size=args.oligomer_size, peak=args.peak, baseline=args.baseline)
        create_tile_layout_viz_from_fasta(args, args.fasta, output_dir, renderer)
        sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="Creates visualizations of FASTA formatted DNA nucleotide data.",
                                     add_help=True)

    parser = argparse.ArgumentParser(prog='Sequenaut.exe')

    parser.add_argument("-f", "--fasta",
                        type=str,
                        help="Path to main FASTA file to process into new visualization.",
                        dest="fasta")
    parser.add_argument("-o", "--output_name",
                        type=str,
                        help="What to name the output folder (not a path). Defaults to name of the fasta file.",
                        dest="output_name")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The layout algorithm.",
                        choices=["triangle_slide", ],
                        dest="layout_type")
    parser.add_argument('--oligomer_size', type=int, default=100, help='Size of the sliding window in nucleotides.',
                        dest="oligomer_size")
    parser.add_argument('--peak', type=float, default=3.0, help='Highest scaling factor for the triangular weight matrix.',
                        dest="peak")
    parser.add_argument('--baseline', type=float, default=1.0, help='Lowest scaling factor for the triangular weight matrix.',
                        dest="baseline")
    args = parser.parse_args()

    if not args.layout_type:
        args.layout_type = "triangle_slide"  # default
    args.output_name = '_'.join([args.output_name, str(args.oligomer_size), str(int(args.peak)), str(int(args.baseline))])

    run_sequenaut(args)
