"""
self.is an optional addon to DDV written in Python that allows you to generate a single image
for an entire genome.  It was necessary to switch platforms and languages because of intrinsic
limitations in the size of image that could be handled by: C#, DirectX, Win2D, GDI+, WIC, SharpDX,
or Direct2D. We tried a lot of options.

self.python file contains basic image handling methods.  It also contains a re-implementation of
Josiah's "Tiled Layout" algorithm which is also in DDVLayoutManager.cs.
"""
import os
import re as regex
import sys
import textwrap

import shutil
from PIL import ImageDraw
import deepzoom




class LayoutLevel:
    def __init__(self, name, modulo, chunk_size=None, padding=0, thickness=1, levels=None):
        self.modulo = modulo
        if chunk_size is not None:
            self.chunk_size = chunk_size
            self._padding = padding
            self.thickness = thickness
        else:
            child = levels[-1]
            self.chunk_size = child.modulo * child.chunk_size
            self._padding = 6 * int(3 ** (len(levels) - 2))  # third level (count=2) should be 6, then 18
            last_parallel = levels[-2]
            self.thickness = last_parallel.modulo * last_parallel.thickness + self.padding

    @property
    def padding(self):
        return self._padding

    @padding.setter
    def padding(self, value):
        original_thickness = self.thickness - self._padding
        self._padding = value
        self.thickness = original_thickness + value


class Contig:
    def __init__(self, name, seq, reset_padding, title_padding, tail_padding, title_index, title_length):
        self.name = name
        self.seq = seq
        self.reset_padding = reset_padding
        self.title_padding = title_padding
        self.tail_padding = tail_padding
        self.nuc_title_start = title_index
        self.nuc_seq_start = title_index + title_length


def multi_line_height(font, multi_line_title, txt):
    sum_line_spacing = ImageDraw.Draw(txt).multiline_textsize(multi_line_title, font)[1]
    descender = font.getsize('y')[1] - font.getsize('A')[1]
    return sum_line_spacing + descender


def pretty_contig_name(contig, title_width, title_lines):
    """Since textwrap.wrap break on whitespace, it's important to make sure there's whitespace
    where there should be.  Contig names don't tend to be pretty."""
    pretty_name = contig.name.replace('_', ' ').replace('|', ' ').replace('chromosome chromosome', 'chromosome')
    pretty_name = regex.sub(r'([^:]*\S):(\S[^:]*)', r'\1: \2', pretty_name)
    pretty_name = regex.sub(r'([^:]*\S):(\S[^:]*)', r'\1: \2', pretty_name)  # don't ask
    if title_width < 20:
        # For small spaces, cram every last bit into the line labels, there's not much room
        pretty_name = pretty_name[:title_width] + '\n' + pretty_name[title_width:title_width * 2]
    else:
        pretty_name = '\n'.join(textwrap.wrap(pretty_name, title_width)[:title_lines])  # approximate width
    return pretty_name


def copytree(src, dst, symlinks=False, ignore=None):
    if not os.path.exists(dst):
        os.makedirs(dst, exist_ok=True)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)


def create_deepzoom_stack(input_image, output_dzi):
    dz_params = {'tile_size': 256,
                 'tile_overlap': 1,
                 'tile_format': "png",
                 'resize_filter': "antialias"}  # cubic bilinear bicubic nearest antialias
    creator = deepzoom.ImageCreator(tile_size=dz_params['tile_size'],
                                    tile_overlap=dz_params['tile_overlap'],
                                    tile_format=dz_params['tile_format'],
                                    resize_filter=dz_params['resize_filter'])
    creator.create(input_image, output_dzi)


def DDV_main(argv):
    from DDVTileLayout import DDVTileLayout
    from ParallelGenomeLayout import ParallelLayout

    folder = "."
    input_file_path = argv[1]
    chromosome_name = os.path.basename(input_file_path[:input_file_path.rfind(".")])  # name between /<path>/ and .png
    image = chromosome_name
    n_arguments = len(argv)

    if n_arguments == 2:  # Shortcut for old visualizations
        output_file = chromosome_name + '.dzi'
        output_dir = os.path.dirname(input_file_path)
        create_deepzoom_stack(input_file_path, os.path.join(output_dir, output_file))
        sys.exit(0)

    if n_arguments == 3:
        raise ValueError("You need to specify an output folder and an image name")

    if n_arguments >= 4:  # Typical use case
        image = argv[3]
        chromosome_name = image  # based on image name, not fasta
        folder = os.path.join(argv[2], chromosome_name)  # place inside a folder with chromosome_name

    if n_arguments > 4:  # Multiple inputs => Parallel genome column layout
        layout = ParallelLayout(n_arguments - 3)
        additional_files = argv[4:]
        layout.process_file(input_file_path, folder, image, additional_files)
    else:  # Typical use case
        layout = DDVTileLayout()
        layout.process_file(input_file_path, folder, image)

    create_deepzoom_stack(os.path.join(folder, image + '.png'), os.path.join(folder, 'GeneratedImages', 'dzc_output.xml'))
    print("Done creating Deep Zoom Structure\nCopying Source File:", input_file_path)
    destination = os.path.join(folder, os.path.basename(input_file_path))
    if not os.path.exists(destination):  # could have been created by ChainParser.py
        shutil.copy(input_file_path, destination)  # copy source file

    # sys.exit(0)


if __name__ == "__main__":
    DDV_main(sys.argv)
