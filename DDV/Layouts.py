import sys


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


def level_layout_factory(modulos, padding, origin):
    # noinspection PyListCreation
    levels = [
        LayoutLevel("XInColumn", modulos[0], 1, padding[0]),  # [0]
        LayoutLevel("LineInColumn", modulos[1], modulos[0], padding[1])  # [1]
    ]
    for i in range(2, len(modulos)):
        levels.append(LayoutLevel("ColumnInRow", modulos[i], padding=padding[i], levels=levels))  # [i]
    return LayoutFrame(origin, levels)


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
