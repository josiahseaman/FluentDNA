import os
import shutil
import re as regex
import textwrap

import deepzoom

from PIL import ImageDraw


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
