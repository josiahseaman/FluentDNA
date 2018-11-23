import sys
from PIL import Image, ImageDraw
from DDV.DDVUtils import multi_line_height


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
            # 6 * int(3 ** (len(levels) - 2))  # third level (count=2) should be 6, then 18
            self._padding = padding or child.padding * 3
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


class LayoutFrame(list):
    """Container class for the origin and LayoutLevels.  There should be one
     LayoutFrame per FASTA source (TileLayout.fasta_sources).
     In every other way this will act like a list containing only the levels."""
    def __init__(self, origin, levels):
        self.origin = origin
        self.levels = levels
        super(LayoutFrame, self).__init__(self.levels)

    @property
    def base_width(self):
        """Shorthand for the column width value that is used often."""
        return self.levels[0].modulo

    def to_json(self):
        return {"origin": list(self.origin), # Origin must be a list, not a tuple
                    "levels": self.levels_json()}

    def levels_json(self):
        json = []
        for level in self.levels:
            json.append({"modulo": level.modulo, "chunk_size": level.chunk_size,
                         "padding": level.padding, "thickness": level.thickness})
        return json

    def relative_position(self, progress):
        """ Readable unoptimized version:
            Maps a nucleotide index to an x,y coordinate based on the rules set in self.levels"""
        xy = [0, 0]
        for i, level in enumerate(self.levels):
            if progress < level.chunk_size:
                return int(xy[0]), int(xy[1])  # somehow a float snuck in here once
            part = i % 2
            coordinate_in_chunk = int(progress // level.chunk_size) % level.modulo
            xy[part] += level.thickness * coordinate_in_chunk
        return [int(xy[0]), int(xy[1])]

    def position_on_screen(self, progress):
        # column padding for various markup = self.levels[2].padding
        xy = self.relative_position(progress)
        return xy[0] + self.origin[0], xy[1] + self.origin[1]


    def write_label(self, contig_name, width, height, font, title_width, upper_left, vertical_label,
                    strand, canvas, horizontal_centering=False, center_vertical=False, chop_text=True,
                    label_color=(50, 50, 50, 255)):
        """write_label() made to nicely draw single line gene labels from annotation
        :param horizontal_centering:
        """
        upper_left = list(upper_left)  # to make it mutable
        shortened = contig_name[-title_width:]  # max length 18.  Last characters are most unique
        txt = Image.new('RGBA', (width, height))#, color=(0,0,0,50))
        txt_canvas = ImageDraw.Draw(txt)
        text_width = txt_canvas.textsize(shortened, font)[0]
        if not chop_text and text_width > width:
            txt = Image.new('RGBA', (text_width, height))  # TODO performance around txt_canvas
            txt_canvas = ImageDraw.Draw(txt)
        if center_vertical or vertical_label:  # Large labels are centered in the column to look nice,
            # rotation indicates strand in big text
            vertically_centered = (height // 2) - multi_line_height(font, shortened, txt)//2
        else:  # Place label at the beginning of gene based on strand
            vertically_centered = height - multi_line_height(font, shortened, txt)  # bottom
            if strand == "+":
                vertically_centered = 0  # top of the box
        txt_canvas.multiline_text((0, max(0, vertically_centered)), shortened, font=font,
                                           fill=label_color)
        if vertical_label:
            rotation_direction = 90 if strand == '-' else -90
            txt = txt.rotate(rotation_direction, expand=True)
            upper_left[1] += -4 if strand == '-' else 4
        if horizontal_centering:
            margin = width - text_width
            upper_left[0] += margin // 2
        canvas.paste(txt, (upper_left[0], upper_left[1]), txt)



def level_layout_factory(modulos, padding, origin):
    # noinspection PyListCreation
    levels = [
        LayoutLevel("XInColumn", modulos[0], 1, padding[0]),  # [0]
        LayoutLevel("LineInColumn", modulos[1], modulos[0], padding[1])  # [1]
    ]
    for i in range(2, len(modulos)):
        levels.append(LayoutLevel("ColumnInRow", modulos[i], padding=padding[i], levels=levels))  # [i]
    return LayoutFrame(tuple(origin), levels)


def parse_custom_layout(custom_layout):
    if custom_layout is not None:
        custom = eval(custom_layout)
        if len(custom) == 2 and hasattr(custom[0], '__iter__') and hasattr(custom[1], '__iter__'):
            modulos, padding = custom
            if all([type(i) == type(7) for i in modulos + padding]) and \
                            len(modulos) == len(padding):
                return modulos, padding
        print('Custom layout must be formatted as two integer lists of euqal length.\n'
                  'For example: --custom_layout="([10,100,100,10,3,999], [0,0,0,3,18,108,200])"', file=sys.stderr)
    return False, False
