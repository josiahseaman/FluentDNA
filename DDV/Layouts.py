

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

    def to_json(self):
        return {"origin": list(self.origin), # Origin must be a list, not a tuple
                    "levels": self.levels_json()}

    def levels_json(self):
        json = []
        for level in self.levels:
            json.append({"modulo": level.modulo, "chunk_size": level.chunk_size,
                         "padding": level.padding, "thickness": level.thickness})
        return json


def level_layout_factory(modulos, padding, origin):
    # noinspection PyListCreation
    levels = [
        LayoutLevel("XInColumn", modulos[0], 1, padding[0]),  # [0]
        LayoutLevel("LineInColumn", modulos[1], modulos[0], padding[1])  # [1]
    ]
    for i in range(2, len(modulos)):
        levels.append(LayoutLevel("ColumnInRow", modulos[i], padding=padding[i], levels=levels))  # [i]
    return LayoutFrame(origin, levels)

