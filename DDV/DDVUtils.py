from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes
import os
import re as regex
import sys
import textwrap
from collections import defaultdict

from PIL import ImageDraw


class keydefaultdict(defaultdict):
    """https://stackoverflow.com/a/2912455/3067894"""
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret


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



def hold_console_for_windows():
    """
    Demonstration of app that knows if it should keep the console open once done.
    Written by Martin Packman.  https://stackoverflow.com/a/2261219/3067894
    """
    try:
        import ctypes
        from ctypes.wintypes import DWORD, HWND

        LPDWORD = ctypes.POINTER(DWORD)

        _kernel32 = ctypes.windll.kernel32
        _user32 = ctypes.windll.user32

        # <http://msdn.microsoft.com/library/ms683180>
        GetCurrentProcessId = ctypes.WINFUNCTYPE(DWORD)(
            ("GetCurrentProcessId", _kernel32))

        # Needs _WIN32_WINNT >= 0x0500
        # <http://msdn.microsoft.com/library/ms683175>
        GetConsoleWindow = ctypes.WINFUNCTYPE(HWND)(("GetConsoleWindow", _kernel32))
        # can return None (NULL)

        # <http://msdn.microsoft.com/library/ms633522>
        GetWindowThreadProcessId = ctypes.WINFUNCTYPE(DWORD, HWND, LPDWORD)(
            ("GetWindowThreadProcessId", _user32),
            ((1, "hWnd"), (2, "lpdwProcessId")))

        # Giving paramflags means the return value gets eaten?

        def owns_console():
            wnd = GetConsoleWindow()
            if wnd is None:
                return False
            return GetCurrentProcessId() == GetWindowThreadProcessId(wnd)


        if owns_console():
            input("Press ENTER key to exit.")
    except BaseException:
        pass  # probably not windows, so it doesn't matter

def viridis_palette():
    """Hard coded copy of Matplotlib's default color palette.  It is
    perceptually uniform and color blind safe."""
    palette = defaultdict(lambda: (255, 0, 0))
    palette[0] = (68, 1, 84)
    palette[1] = (68, 2, 85)
    palette[2] = (68, 3, 87)
    palette[3] = (69, 5, 88)
    palette[4] = (69, 6, 90)
    palette[5] = (69, 8, 91)
    palette[6] = (70, 9, 92)
    palette[7] = (70, 11, 94)
    palette[8] = (70, 12, 95)
    palette[9] = (70, 14, 97)
    palette[10] = (71, 15, 98)
    palette[11] = (71, 17, 99)
    palette[12] = (71, 18, 101)
    palette[13] = (71, 20, 102)
    palette[14] = (71, 21, 103)
    palette[15] = (71, 22, 105)
    palette[16] = (71, 24, 106)
    palette[17] = (72, 25, 107)
    palette[18] = (72, 26, 108)
    palette[19] = (72, 28, 110)
    palette[20] = (72, 29, 111)
    palette[21] = (72, 30, 112)
    palette[22] = (72, 32, 113)
    palette[23] = (72, 33, 114)
    palette[24] = (72, 34, 115)
    palette[25] = (72, 35, 116)
    palette[26] = (71, 37, 117)
    palette[27] = (71, 38, 118)
    palette[28] = (71, 39, 119)
    palette[29] = (71, 40, 120)
    palette[30] = (71, 42, 121)
    palette[31] = (71, 43, 122)
    palette[32] = (71, 44, 123)
    palette[33] = (70, 45, 124)
    palette[34] = (70, 47, 124)
    palette[35] = (70, 48, 125)
    palette[36] = (70, 49, 126)
    palette[37] = (69, 50, 127)
    palette[38] = (69, 52, 127)
    palette[39] = (69, 53, 128)
    palette[40] = (69, 54, 129)
    palette[41] = (68, 55, 129)
    palette[42] = (68, 57, 130)
    palette[43] = (67, 58, 131)
    palette[44] = (67, 59, 131)
    palette[45] = (67, 60, 132)
    palette[46] = (66, 61, 132)
    palette[47] = (66, 62, 133)
    palette[48] = (66, 64, 133)
    palette[49] = (65, 65, 134)
    palette[50] = (65, 66, 134)
    palette[51] = (64, 67, 135)
    palette[52] = (64, 68, 135)
    palette[53] = (63, 69, 135)
    palette[54] = (63, 71, 136)
    palette[55] = (62, 72, 136)
    palette[56] = (62, 73, 137)
    palette[57] = (61, 74, 137)
    palette[58] = (61, 75, 137)
    palette[59] = (61, 76, 137)
    palette[60] = (60, 77, 138)
    palette[61] = (60, 78, 138)
    palette[62] = (59, 80, 138)
    palette[63] = (59, 81, 138)
    palette[64] = (58, 82, 139)
    palette[65] = (58, 83, 139)
    palette[66] = (57, 84, 139)
    palette[67] = (57, 85, 139)
    palette[68] = (56, 86, 139)
    palette[69] = (56, 87, 140)
    palette[70] = (55, 88, 140)
    palette[71] = (55, 89, 140)
    palette[72] = (54, 90, 140)
    palette[73] = (54, 91, 140)
    palette[74] = (53, 92, 140)
    palette[75] = (53, 93, 140)
    palette[76] = (52, 94, 141)
    palette[77] = (52, 95, 141)
    palette[78] = (51, 96, 141)
    palette[79] = (51, 97, 141)
    palette[80] = (50, 98, 141)
    palette[81] = (50, 99, 141)
    palette[82] = (49, 100, 141)
    palette[83] = (49, 101, 141)
    palette[84] = (49, 102, 141)
    palette[85] = (48, 103, 141)
    palette[86] = (48, 104, 141)
    palette[87] = (47, 105, 141)
    palette[88] = (47, 106, 141)
    palette[89] = (46, 107, 142)
    palette[90] = (46, 108, 142)
    palette[91] = (46, 109, 142)
    palette[92] = (45, 110, 142)
    palette[93] = (45, 111, 142)
    palette[94] = (44, 112, 142)
    palette[95] = (44, 113, 142)
    palette[96] = (44, 114, 142)
    palette[97] = (43, 115, 142)
    palette[98] = (43, 116, 142)
    palette[99] = (42, 117, 142)
    palette[100] = (42, 118, 142)
    palette[101] = (42, 119, 142)
    palette[102] = (41, 120, 142)
    palette[103] = (41, 121, 142)
    palette[104] = (40, 122, 142)
    palette[105] = (40, 122, 142)
    palette[106] = (40, 123, 142)
    palette[107] = (39, 124, 142)
    palette[108] = (39, 125, 142)
    palette[109] = (39, 126, 142)
    palette[110] = (38, 127, 142)
    palette[111] = (38, 128, 142)
    palette[112] = (38, 129, 142)
    palette[113] = (37, 130, 142)
    palette[114] = (37, 131, 141)
    palette[115] = (36, 132, 141)
    palette[116] = (36, 133, 141)
    palette[117] = (36, 134, 141)
    palette[118] = (35, 135, 141)
    palette[119] = (35, 136, 141)
    palette[120] = (35, 137, 141)
    palette[121] = (34, 137, 141)
    palette[122] = (34, 138, 141)
    palette[123] = (34, 139, 141)
    palette[124] = (33, 140, 141)
    palette[125] = (33, 141, 140)
    palette[126] = (33, 142, 140)
    palette[127] = (32, 143, 140)
    palette[128] = (32, 144, 140)
    palette[129] = (32, 145, 140)
    palette[130] = (31, 146, 140)
    palette[131] = (31, 147, 139)
    palette[132] = (31, 148, 139)
    palette[133] = (31, 149, 139)
    palette[134] = (31, 150, 139)
    palette[135] = (30, 151, 138)
    palette[136] = (30, 152, 138)
    palette[137] = (30, 153, 138)
    palette[138] = (30, 153, 138)
    palette[139] = (30, 154, 137)
    palette[140] = (30, 155, 137)
    palette[141] = (30, 156, 137)
    palette[142] = (30, 157, 136)
    palette[143] = (30, 158, 136)
    palette[144] = (30, 159, 136)
    palette[145] = (30, 160, 135)
    palette[146] = (31, 161, 135)
    palette[147] = (31, 162, 134)
    palette[148] = (31, 163, 134)
    palette[149] = (32, 164, 133)
    palette[150] = (32, 165, 133)
    palette[151] = (33, 166, 133)
    palette[152] = (33, 167, 132)
    palette[153] = (34, 167, 132)
    palette[154] = (35, 168, 131)
    palette[155] = (35, 169, 130)
    palette[156] = (36, 170, 130)
    palette[157] = (37, 171, 129)
    palette[158] = (38, 172, 129)
    palette[159] = (39, 173, 128)
    palette[160] = (40, 174, 127)
    palette[161] = (41, 175, 127)
    palette[162] = (42, 176, 126)
    palette[163] = (43, 177, 125)
    palette[164] = (44, 177, 125)
    palette[165] = (46, 178, 124)
    palette[166] = (47, 179, 123)
    palette[167] = (48, 180, 122)
    palette[168] = (50, 181, 122)
    palette[169] = (51, 182, 121)
    palette[170] = (53, 183, 120)
    palette[171] = (54, 184, 119)
    palette[172] = (56, 185, 118)
    palette[173] = (57, 185, 118)
    palette[174] = (59, 186, 117)
    palette[175] = (61, 187, 116)
    palette[176] = (62, 188, 115)
    palette[177] = (64, 189, 114)
    palette[178] = (66, 190, 113)
    palette[179] = (68, 190, 112)
    palette[180] = (69, 191, 111)
    palette[181] = (71, 192, 110)
    palette[182] = (73, 193, 109)
    palette[183] = (75, 194, 108)
    palette[184] = (77, 194, 107)
    palette[185] = (79, 195, 105)
    palette[186] = (81, 196, 104)
    palette[187] = (83, 197, 103)
    palette[188] = (85, 198, 102)
    palette[189] = (87, 198, 101)
    palette[190] = (89, 199, 100)
    palette[191] = (91, 200, 98)
    palette[192] = (94, 201, 97)
    palette[193] = (96, 201, 96)
    palette[194] = (98, 202, 95)
    palette[195] = (100, 203, 93)
    palette[196] = (103, 204, 92)
    palette[197] = (105, 204, 91)
    palette[198] = (107, 205, 89)
    palette[199] = (109, 206, 88)
    palette[200] = (112, 206, 86)
    palette[201] = (114, 207, 85)
    palette[202] = (116, 208, 84)
    palette[203] = (119, 208, 82)
    palette[204] = (121, 209, 81)
    palette[205] = (124, 210, 79)
    palette[206] = (126, 210, 78)
    palette[207] = (129, 211, 76)
    palette[208] = (131, 211, 75)
    palette[209] = (134, 212, 73)
    palette[210] = (136, 213, 71)
    palette[211] = (139, 213, 70)
    palette[212] = (141, 214, 68)
    palette[213] = (144, 214, 67)
    palette[214] = (146, 215, 65)
    palette[215] = (149, 215, 63)
    palette[216] = (151, 216, 62)
    palette[217] = (154, 216, 60)
    palette[218] = (157, 217, 58)
    palette[219] = (159, 217, 56)
    palette[220] = (162, 218, 55)
    palette[221] = (165, 218, 53)
    palette[222] = (167, 219, 51)
    palette[223] = (170, 219, 50)
    palette[224] = (173, 220, 48)
    palette[225] = (175, 220, 46)
    palette[226] = (178, 221, 44)
    palette[227] = (181, 221, 43)
    palette[228] = (183, 221, 41)
    palette[229] = (186, 222, 39)
    palette[230] = (189, 222, 38)
    palette[231] = (191, 223, 36)
    palette[232] = (194, 223, 34)
    palette[233] = (197, 223, 33)
    palette[234] = (199, 224, 31)
    palette[235] = (202, 224, 30)
    palette[236] = (205, 224, 29)
    palette[237] = (207, 225, 28)
    palette[238] = (210, 225, 27)
    palette[239] = (212, 225, 26)
    palette[240] = (215, 226, 25)
    palette[241] = (218, 226, 24)
    palette[242] = (220, 226, 24)
    palette[243] = (223, 227, 24)
    palette[244] = (225, 227, 24)
    palette[245] = (228, 227, 24)
    palette[246] = (231, 228, 25)
    palette[247] = (233, 228, 25)
    palette[248] = (236, 228, 26)
    palette[249] = (238, 229, 27)
    palette[250] = (241, 229, 28)
    palette[251] = (243, 229, 30)
    palette[252] = (246, 230, 31)
    palette[253] = (248, 230, 33)
    palette[254] = (250, 230, 34)
    palette[255] = (253, 231, 36)
    return palette