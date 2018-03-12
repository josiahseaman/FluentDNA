from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import argparse
import math
import sys
from collections import defaultdict

from DDV.TileLayout import TileLayout
from DDV.fluentdna import create_tile_layout_viz_from_fasta
from DDV.DDVUtils import make_output_dir_with_suffix, base_directories


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


def get_line(start, end):
    """Bresenham's Line Algorithm
    Produces a list of tuples from start and end
    Copied from http://www.roguebasin.com/index.php?title=Bresenham%27s_Line_Algorithm#Python
    """
    # Setup initial conditions
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1

    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1

    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1

    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()
    return points


class Sequenaut(TileLayout):  # TODO: make an abstract class parent
    def __init__(self, layout='triangle_slide', oligomer_size=100, peak=10.0, baseline=1.0, log_scale=True):
        super(Sequenaut, self).__init__()
        self.oligomer_size = oligomer_size
        self.layout = layout
        self.peak = peak
        self.baseline = baseline
        midpoint = self.oligomer_size // 2
        self.weight_matrix = list(linspace(baseline, self.peak, midpoint)) + list(linspace(self.peak, baseline, self.oligomer_size - midpoint))
        self.image_length = self.max_dimensions(oligomer_size, verbose=True)[0]
        self.hits = [[0] * self.image_length for i in range(self.image_length)]
        self.log_scale = log_scale

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


    def connected_sequenaut(self, seq):
        prev_coord = self.coordinate(seq[: self.oligomer_size], self.weight_matrix)
        for i in range(len(seq) - self.oligomer_size + 1):
            olig = seq[i: i + self.oligomer_size]
            point = self.coordinate(olig, self.weight_matrix)
            line = get_line(prev_coord, point)
            for x, y in line:
                self.hits[y][x] += 1
            prev_coord = point


    def draw_nucleotides(self):
        for contig in self.contigs:
            if self.layout == 'connected':
                self.connected_sequenaut(contig.seq)
            else:
                self.weighted_sequenaut(contig.seq)
        leader = max([max(column) for column in self.hits])
        if self.log_scale:
            leader = math.log(leader)
        try:
            middle = self.max_dimensions(self.oligomer_size)[0] // 2
            print("Leader", leader, int(math.log(self.hits[middle][middle]) / leader * 235 + 20))
        except: pass  # just talking

        for y in range(len(self.hits)):
            for x in range(len(self.hits[y])):
                if self.hits[y][x]:
                    val = math.log(self.hits[y][x]) if self.log_scale else self.hits[y][x]
                    grey = 255 - int(val / leader * 235 + 20)
                    self.pixels[x, y] = (grey, grey, grey)

    def draw_titles(self):
        pass  # There are no titles in Sequenaut


    def process_file(self, input_file_path, output_folder, output_file_name):
        # use the parent process_file() but override methods in child
        super(Sequenaut, self).process_file(input_file_path, output_folder, output_file_name)
        # TODO: support --no_webpage
        self.generate_html(output_folder, output_file_name)


def run_sequenaut(args):
    SERVER_HOME, base_path = base_directories(args)
    # TODO: allow batch of tiling layout by chromosome
    output_dir = make_output_dir_with_suffix(base_path, '')
    renderer = Sequenaut(layout=args.layout, oligomer_size=args.oligomer_size, peak=args.peak, baseline=args.baseline, log_scale=not args.linear_scale)
    create_tile_layout_viz_from_fasta(args, args.fasta, output_dir, args.output_name, renderer)
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
                        choices=["triangle_slide", "connected"],)
    parser.add_argument('--linear_scale', action='store_true', help='Use linear color scaling (defaults to logarithmic).',)
    parser.add_argument('--oligomer_size', type=int, default=100, help='Size of the sliding window in nucleotides.',)
    parser.add_argument('--peak', type=float, default=3.0, help='Highest scaling factor for the triangular weight matrix.',)
    parser.add_argument('--baseline', type=float, default=1.0, help='Lowest scaling factor for the triangular weight matrix.',)
    args = parser.parse_args()

    if not args.layout:
        args.layout = "triangle_slide"  # default
    args.output_name = '_'.join([args.output_name, str(args.oligomer_size), str(int(args.peak)), str(int(args.baseline))])

    run_sequenaut(args)
