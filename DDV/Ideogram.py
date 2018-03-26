"""Plot a DNA sequence along a murray space filling curve.

This FluentDNA Python class is based on Yan Wong's code (unpublished) which is in turn
based on: algol version from https://alb.host.cs.st-andrews.ac.uk/cole/code.html
(see also http://www.aleph.se/andart/archives/2013/10/murray_polygons.html)

The murray polygon is stretched by 3 times in the X and 5 times in the Y direction
to improve the locality properties (ideally this should be 1:sqrt(3) = 0.577 rather
than 3:5 = 0.6) see https://arxiv.org/pdf/0806.4787.pdf

try e.g. python3 Ideogram.py -x 3 3 3 -y 3 3 3
"""


from DDV.TileLayout import TileLayout, hex_to_rgb
from PIL import Image, ImageDraw
import os
import numpy as np

class Ideogram(TileLayout):
    def __init__(self, x_radices, y_radices):
        self.x_radices = x_radices
        self.y_radices = y_radices
        super(Ideogram, self).__init__()


    def draw_nucleotides(self):
        ndim_x = len(self.x_radices)
        ndim_y = len(self.y_radices)
        max_dim = max(ndim_x, ndim_y)

        radices = np.ones((max_dim, 2), dtype=np.int)
        digits = np.zeros((max_dim * 2), dtype=np.int)
        parities = np.ones((max_dim + 1, 2), dtype=np.int)
        prevprev_pos = np.zeros((2,), dtype=np.int)  # stores y,x
        prev_pos = np.zeros((2,), dtype=np.int)  # stores y,x
        curr_pos = np.zeros((2,), dtype=np.int)  # stores y,x
        radices[0:ndim_x, 0] = self.x_radices
        radices[0:ndim_y, 1] = self.y_radices

        no_pts = np.prod(radices)
        points_visited = set()
        n = np.prod(radices, 0)  # get nx = n[0] and ny = n[1]
        radices.shape = np.prod(radices.shape)  # flatten

        points_file_name = os.path.join(self.final_output_location, "test_ideogram_points.txt")
        points_file = None # open(points_file_name, 'w')
        if points_file:
            print("Saving locations in {}".format(points_file_name))
        contig = self.contigs[0]
        seq_iter = iter(contig.seq)
        x_scale, y_scale = 1, 1

        if x_scale == 1 and y_scale == 1:
            self.draw_loop_optimized( curr_pos, digits, no_pts, parities, points_visited, prev_pos, radices,
                                      seq_iter)
        else:
            self.draw_loop_any_scale(curr_pos, digits, no_pts, parities, points_file, points_visited,
                                     prev_pos, prevprev_pos, radices, seq_iter, x_scale, y_scale)


    def draw_loop_optimized(self, curr_pos, digits, no_pts, parities, points_visited, prev_pos, radices,
                            seq_iter):
        for pts in range(no_pts - 1):
            place = increment(digits, radices, 0)
            parities[0:(place // 2 + 1), place % 2] *= -1
            place += 1
            prev_pos[:] = curr_pos[:]
            # assume we move 3 up and 5 across
            x = int(prev_pos[1] + 2)
            y = int(prev_pos[0] + 2)
            curr_pos[place % 2] += parities[place // 2, place % 2]
            assert (x, y) not in points_visited
            points_visited.add((x, y))
            try:
                self.draw_pixel(next(seq_iter), x, y)
            except StopIteration:
                break  # reached end of sequence


    def draw_loop_any_scale(self, curr_pos, digits, no_pts, parities, points_file, points_visited, prev_pos,
                            prevprev_pos, radices, seq_iter, x_scale, y_scale):
        for pts in range(no_pts - 1):
            place = increment(digits, radices, 0)
            parities[0:(place // 2 + 1), place % 2] *= -1
            place += 1
            prevprev_pos[:] = prev_pos[:]
            prev_pos[:] = curr_pos[:]
            curr_pos[place % 2] += parities[place // 2, place % 2]
            # assume we move 3 up and 5 across
            x = int(prev_pos[1] * x_scale + 2)
            y = int(prev_pos[0] * y_scale + 2)
            if points_file:
                print("{} {}".format(x, y), file=points_file)
            diff = curr_pos - prev_pos
            prev_diff = prev_pos - prevprev_pos
            assert (abs(sum(diff)) == 1)
            assert (x, y) not in points_visited
            points_visited.add((x, y))
            try:
                self.paint_turns(seq_iter, x, y, diff, prev_diff,
                                 prev_pos, prevprev_pos, x_scale, y_scale)
            except IndexError:
                print(x, y, "out of range")
            except StopIteration:
                break  # reached end of sequence


    def paint_turns(self, seq_iter, x, y, diff, prev_diff, prev_pos, prevprev_pos, x_scale, y_scale):
        # right-hand rotation at corner when corner==1, left-hand rotation when corner==1, or no turn (corner == 0)
        turn = prev_diff[0] * diff[1] - prev_diff[1] * diff[0]
        if turn == 0 or (x_scale == 1 and y_scale==1):
            self.draw_pixel(next(seq_iter), x, y)
        if diff[1]:
            # x is changing
            for scale_step in range(1, x_scale):
                self.draw_pixel(next(seq_iter), x + scale_step * int(diff[1]), y)
        elif diff[0]:
            # y is changing
            # NB: underlines will sometimes overwrite previous ones
            x_nudge = int(prevprev_pos[1] - prev_pos[1])
            for scale_step in range(1, y_scale):
                self.draw_pixel(next(seq_iter), x + x_nudge, y + scale_step * int(diff[0]))

    def max_dimensions(self, image_length):
        dim = int(np.sqrt(image_length * 2))  # ideogram has low density and mostly square
        return dim, dim


    def prepare_image(self, image_length):
        ui_grey = hex_to_rgb('EEF3FA')  # The contrast isn't quite so bad as white
        width, height = self.max_dimensions(image_length)
        print("Image dimensions are", width, "x", height, "pixels")
        self.image = Image.new('RGB', (width, height), ui_grey)
        self.draw = ImageDraw.Draw(self.image)
        self.pixels = self.image.load()

def increment(digits, radices, place):
    if digits[place] < (radices[place] - 1):
        digits[place] += 1
        return place
    else:
        digits[place] = 0
        return increment(digits,radices,place + 1)


if __name__ == "__main__":
    layout = Ideogram([3,3,3,3,3,3,3], [5,5,3,3,3,3])
    layout.process_file("example_data/hg38_chr19_sample.fa", 'www-data/dnadata/test ideogram', 'ideogram')
