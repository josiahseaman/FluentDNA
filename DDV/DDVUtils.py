from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes
import os
import re as regex
import sys
import textwrap

from PIL import ImageDraw


class LayoutLevel(object):
    def __init__(self, name, modulo, chunk_size=None, padding=None, thickness=1, levels=None):
        self.modulo = modulo
        if chunk_size is not None:
            self.chunk_size = chunk_size
            self._padding = padding
            self.thickness = thickness
        else:
            child = levels[-1]
            self.chunk_size = child.modulo * child.chunk_size
            self._padding = padding or child.padding * 3  # 6 * int(3 ** (len(levels) - 2))  # third level (count=2) should be 6, then 18
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



def multi_line_height(font, multi_line_title, txt):
    sum_line_spacing = ImageDraw.Draw(txt).multiline_textsize(multi_line_title, font)[1]
    descender = font.getsize('y')[1] - font.getsize('A')[1]
    return sum_line_spacing + descender


def pretty_contig_name(contig_name, title_width, title_lines):
    """Since textwrap.wrap break on whitespace, it's important to make sure there's whitespace
    where there should be.  Contig names don't tend to be pretty."""
    pretty_name = contig_name.replace('_', ' ').replace('|', ' ').replace('chromosome chromosome', 'chromosome')
    pretty_name = regex.sub(r'([^:]*\S):(\S[^:]*)', r'\1: \2', pretty_name)
    pretty_name = regex.sub(r'([^:]*\S):(\S[^:]*)', r'\1: \2', pretty_name)  # don't ask
    if title_width < 20 and len(pretty_name) > title_width * 1.5:  # this is a suboptimal special case to try and
        # cram more characters onto the two lines of the smallest contig titles when there's not enough space
        # For small spaces, cram every last bit into the line labels, there's not much room
        pretty_name = pretty_name[:title_width] + '\n' + pretty_name[title_width:title_width * 2]
    else:  # this is the only case that correctly bottom justifies one line titles
        pretty_name = '\n'.join(textwrap.wrap(pretty_name, title_width)[:title_lines])  # approximate width
    return pretty_name


def create_deepzoom_stack(input_image, output_dzi):
    import DDV.deepzoom as deepzoom
    creator = deepzoom.ImageCreator(tile_size=256,
                                    tile_overlap=1,
                                    tile_format="png",
                                    resize_filter="antialias")# cubic bilinear bicubic nearest antialias
    creator.create(input_image, output_dzi)


def make_output_dir_with_suffix(base_path, suffix):
    from os import errno
    output_dir = base_path + suffix
    print("Creating Chromosome Output Directory...", os.path.basename(output_dir))
    try:
        os.makedirs(output_dir)
    except OSError as e:  # exist_ok=True
        if e.errno != errno.EEXIST:
            raise
    return output_dir


def base_directories(args):
    if getattr(sys, 'frozen', False):
        BASE_DIR = os.path.dirname(sys.executable)
    else:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    SERVER_HOME = os.path.join(BASE_DIR, 'www-data', 'dnadata')
    base_path = os.path.join(SERVER_HOME, args.output_name) if args.output_name else SERVER_HOME
    return SERVER_HOME, base_path

