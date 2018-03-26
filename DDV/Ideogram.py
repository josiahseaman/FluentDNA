"""Plot a DNA sequence along a murray space filling curve.

This FluentDNA Python class is based on Yan Wong's code (unpublished) which is in turn
based on: algol version from https://alb.host.cs.st-andrews.ac.uk/cole/code.html
(see also http://www.aleph.se/andart/archives/2013/10/murray_polygons.html)

The murray polygon is stretched by 3 times in the X and 5 times in the Y direction
to improve the locality properties (ideally this should be 1:sqrt(3) = 0.577 rather
than 3:5 = 0.6) see https://arxiv.org/pdf/0806.4787.pdf

try e.g. python3 Ideogram.py -x 3 3 3 -y 3 3 3
"""


from DDV.TileLayout import TileLayout
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
        n = np.prod(radices, 0)  # get nx = n[0] and ny = n[1]
        radices.shape = np.prod(radices.shape)  # flatten

        points_file_name = "data/test_ideogram_points.txt"
        points_file = open(points_file_name, 'w')
        if points_file:
            print("Saving locations in {}".format(points_file_name))
        contig = self.contigs[0]
        seq_iter = iter(contig.seq)

        for pts in range(no_pts - 1):
            place = increment(digits, radices, 0)
            parities[0:(place // 2 + 1), place % 2] *= -1
            place += 1
            prevprev_pos[:] = prev_pos[:]
            prev_pos[:] = curr_pos[:]
            curr_pos[place % 2] += parities[place // 2, place % 2]
            # assume we move 3 up and 5 across
            x = int(prev_pos[1] * 3 + 1)
            y = int(prev_pos[0] * 5 + 2)
            if points_file:
                print("{} {}".format(x, y), file=points_file)
            diff = curr_pos - prev_pos
            prev_diff = prev_pos - prevprev_pos
            assert (abs(sum(diff)) == 1)
            self.paint_turns(seq_iter, x, y, diff, prev_diff, prev_pos, prevprev_pos)

    def paint_turns(self, seq_iter, x, y, diff, prev_diff, prev_pos, prevprev_pos):
        # right-hand rotation at corner when corner==1, left-hand rotation when corner==1, or no turn (corner == 0)
        turn = prev_diff[0] * diff[1] - prev_diff[1] * diff[0]
        if turn == 0:
            self.draw_pixel(next(seq_iter), x, y)
        if diff[1]:
            # x is changing
            self.draw_pixel(next(seq_iter), x + int(diff[1]), y)
            self.draw_pixel(next(seq_iter), x + 2 * int(diff[1]), y)
        else:
            # y is changing
            # NB: underlines will sometimes overwrite previous ones
            x_nudge = int(prevprev_pos[1] - prev_pos[1])
            self.draw_pixel(next(seq_iter), x + x_nudge, y + int(diff[0]))
            self.draw_pixel(next(seq_iter), x + x_nudge, y + 2 * int(diff[0]))
            self.draw_pixel(next(seq_iter), x + x_nudge, y + 3 * int(diff[0]))
            self.draw_pixel(next(seq_iter), x + x_nudge, y + 4 * int(diff[0]))

    def max_dimensions(self, image_length):
        dim = int(np.sqrt(image_length * 11))  # ideogram has low density and mostly square
        return dim, dim


def increment(digits, radices, place):
    if digits[place] < (radices[place] - 1):
        digits[place] += 1
        return place
    else:
        digits[place] = 0
        return increment(digits,radices,place + 1)


if __name__ == "__main__":
    layout = Ideogram([3,3,3], [3,3,3])
    layout.process_file("example_data/phiX.fa", 'www-data/dnadata/test ideogram', 'ideogram')
