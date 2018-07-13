from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import math
import os
import traceback
from collections import defaultdict
from datetime import datetime
from PIL import Image, ImageDraw, ImageFont

from DNASkittleUtils.Contigs import read_contigs, Contig
from DNASkittleUtils.DDVUtils import copytree
from DDV.DDVUtils import LayoutLevel, multi_line_height, pretty_contig_name, viridis_palette
from DDV import gap_char

small_title_bp = 10000
font_name = "Arial.ttf"
try:
    ImageFont.truetype(font_name, 10)
except IOError:
    try:
        font_name = font_name.lower()  # windows and mac are both case sensitive in opposite directions
        ImageFont.truetype(font_name, 10)
    except IOError:
        font_name = None



def hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2 ,4))


def level_layout_factory(modulos, padding=None):
    if padding is None:
        padding = [0, 0, 6, 6 * 3, 6 * (3 ** 2), 6 * (3 ** 3), 6 * (3 ** 4)]
    # noinspection PyListCreation
    levels = [
        LayoutLevel("XInColumn", modulos[0], 1, padding[0]),  # [0]
        LayoutLevel("LineInColumn", modulos[1], modulos[0], padding[1])  # [1]
    ]
    for i in range(2, len(modulos)):
        levels.append(LayoutLevel("ColumnInRow", modulos[i], padding=padding[i], levels=levels))  # [i]
    # levels.append(LayoutLevel("RowInTile", modulos[3], levels=levels))  # [3]
    # levels.append(LayoutLevel("TileColumn", modulos[4], levels=levels))  # [4]
    # levels.append(LayoutLevel("TileRow", modulos[5], levels=levels))  # [5]
    # levels.append(LayoutLevel("PageColumn", modulos[6], levels=levels))  # [6]
    return levels


class TileLayout(object):

    def __init__(self, use_fat_headers=False, use_titles=True, sort_contigs=False,
                 low_contrast=False, base_width=None):
        # use_fat_headers: For large chromosomes in multipart files, do you change the layout to allow for titles that
        # are outside of the nucleotide coordinate grid?
        self.use_titles = use_titles
        self.use_fat_headers = use_fat_headers  # Can only be changed in code.
        self.skip_small_titles = False
        self.using_spectrum = False
        self.sort_contigs = sort_contigs
        self.low_contrast = low_contrast
        self.base_width = base_width if base_width is not None else 100
        self.title_skip_padding = self.base_width  # skip one line. USER: Change this

        # precomputing fonts turns out to be a big performance gain
        sizes = [9, 38, 380, 380 * 2]
        if font_name is not None:
            self.fonts = {size: ImageFont.truetype(font_name, size) for size in sizes}
        else:
            self.fonts = {size: ImageFont.load_default() for size in sizes}
        self.final_output_location = None
        self.image = None
        self.draw = None
        self.pixels = None
        self.contigs = []
        self.image_length = 0
        #Natural, color blind safe Colors
        self.palette = defaultdict(lambda: (255, 0, 0))  # default red will stand out

        #### Rasmol Protein colors
        self.palette['D'] = hex_to_rgb('EA3535')
        self.palette['E'] = hex_to_rgb('EA3535')
        self.palette['F'] = hex_to_rgb('4B4BB5')
        self.palette['H'] = hex_to_rgb('9595D9')
        self.palette['I'] = hex_to_rgb('2F932F')
        self.palette['K'] = hex_to_rgb('3C76FF')
        self.palette['L'] = hex_to_rgb('2F932F')
        self.palette['M'] = hex_to_rgb('ECEC41')
        self.palette['N'] = hex_to_rgb('3BE4E4')
        self.palette['P'] = hex_to_rgb('E25826')
        self.palette['Q'] = hex_to_rgb('3BE4E4')
        self.palette['R'] = hex_to_rgb('3C76FF')
        self.palette['S'] = hex_to_rgb('FBAC34')
        self.palette['V'] = hex_to_rgb('2F932F')
        self.palette['W'] = hex_to_rgb('BF72BF')
        self.palette['X'] = hex_to_rgb('FF6100')
        self.palette['Y'] = hex_to_rgb('4B4BB5')

        self.palette['N'] = (61, 61, 61)  # charcoal grey
        self.palette[gap_char] = (247, 247, 247)  # almost white
        self.palette['.'] = self.palette[gap_char]  # other gap characters

        self.activate_high_contrast_colors()
        if self.low_contrast:
            self.activate_natural_colors()

        # Used in translocations, not amino acids:  B J O U Z
        self.palette['-'] = self.palette[gap_char]
        self.palette['J'] = hex_to_rgb('#E5F3FF')  #E5F3FF blue
        self.palette['B'] = hex_to_rgb('#FFF0EF')  #EAFFE5 green
        self.palette['O'] = hex_to_rgb('#FFEEED')  #FFE7E5 red
        self.palette['Z'] = hex_to_rgb('#F9EDFF')  #F8E5FF pink
        self.palette['U'] = hex_to_rgb('#FFF3E5')  #FFF3E5 orange

        modulos = [self.base_width, self.base_width * 10, 100, 10, 3, 4, 999]
        padding = [0, 0, 6, 6*3, 6*(3**2), 6*(3**3), 6*(3**4)]
        self.levels = level_layout_factory(modulos, padding=padding)

        self.tile_label_size = self.levels[3].chunk_size * 2
        self.origin = [self.levels[2].padding, self.levels[2].padding]
        if self.use_fat_headers:
            self.enable_fat_headers()



    def activate_high_contrast_colors(self):
        # Original DDV Colors
        self.palette['A'] = (255, 0, 0)
        self.palette['G'] = (0, 255, 0)
        self.palette['T'] = (250, 240, 114)
        self.palette['C'] = (0, 0, 255)

    def activate_natural_colors(self):
        # -----Nucleotide Colors! Paletton Quadrapole colors------
        # self.palette['A'] = hex_to_rgb('C35653')  # Red
        # self.palette['T'] = hex_to_rgb('D4A16A')  # Yellow
        # self.palette['G'] = hex_to_rgb('55AA55')  # Green
        # self.palette['C'] = hex_to_rgb('457585')  # Blue
        # -----Nucleotide Colors! Paletton darks ------
        # self.palette['A'] = hex_to_rgb('B94A24')  # Red
        # self.palette['T'] = hex_to_rgb('B98124')  # Yellow
        # self.palette['G'] = hex_to_rgb('19814F')  # Green
        # self.palette['C'] = hex_to_rgb('20467A')  # Blue
        # # -----Nucleotide Colors! Paletton Pastel------
        # self.palette['A'] = hex_to_rgb('EC8D6C')  # Red
        # self.palette['T'] = hex_to_rgb('ECBC6C')  # Yellow
        # self.palette['G'] = hex_to_rgb('4CA47A')  # Green
        # self.palette['C'] = hex_to_rgb('4F6F9B')  # Blue
        # # -----Nucleotide Colors! Paletton Stark ------
        # self.palette['A'] = hex_to_rgb('FF4100')  # Red
        # self.palette['T'] = hex_to_rgb('FF9F00')  # Yellow
        # self.palette['G'] = hex_to_rgb('00C566')  # Green
        # self.palette['C'] = hex_to_rgb('0B56BE')  # Blue
        # -----Manually Adjusted Colors from Paletton plus contrast------
        self.palette['A'] = hex_to_rgb('D4403C')  # Red
        self.palette['T'] = hex_to_rgb('E2AE5B')  # Yellow
        self.palette['G'] = hex_to_rgb('3FB93F')  # Green
        self.palette['C'] = hex_to_rgb('2D6C85')  # Blue

    def enable_fat_headers(self):
        if self.use_titles:
            print("Using Fat Headers!")
            self.use_fat_headers = True
            self.levels = self.levels[:6]
            self.levels[5].padding += self.levels[3].thickness  # one full row for a chromosome title
            self.levels.append(LayoutLevel("PageColumn", 999, levels=self.levels))  # [6]
            self.origin[1] += self.levels[5].padding  # padding comes before, not after
            self.tile_label_size = 0  # Fat_headers are not part of the coordinate space

    def process_file(self, input_file_path, output_folder, output_file_name):
        start_time = datetime.now()
        self.image_length = self.read_contigs_and_calc_padding(input_file_path)
        print("Read contigs :", datetime.now() - start_time)
        self.prepare_image(self.image_length)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n')
            traceback.print_exc()
        try:
            if len(self.contigs) > 1 and self.use_titles:
                print("Drawing %i titles" % sum(len(x.seq) > small_title_bp for x in self.contigs))
                self.draw_titles()
                print("Drew Titles:", datetime.now() - start_time)
        except BaseException as e:
            print('Encountered exception while drawing titles:', '\n')
            traceback.print_exc()
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)

    def draw_nucleotides(self):
        total_progress = 0
        # Layout contigs one at a time
        for contig in self.contigs:
            total_progress += contig.reset_padding + contig.title_padding
            seq_length = len(contig.seq)
            line_width = self.levels[0].modulo
            for cx in range(0, seq_length, line_width):
                x, y = self.position_on_screen(total_progress)
                remaining = min(line_width, seq_length - cx)
                total_progress += remaining
                #try:
                for i in range(remaining):
                    nuc = contig.seq[cx + i]
                    # if nuc != gap_char:
                    self.draw_pixel(nuc, x + i, y)
                #except IndexError:
                #    print("Cursor fell off the image at", x,y)
                if cx % 100000 == 0:
                    print('\r', str(total_progress / self.image_length * 100)[:6], '% done:', contig.name,
                          end="")  # pseudo progress bar
            total_progress += contig.tail_padding  # add trailing white space after the contig sequence body
        print('')

    def calc_all_padding(self):
        total_progress = 0  # pointer in image
        seq_start = 0  # pointer in text
        multipart_file = len(self.contigs) > 1

        # if len(self.levels) >= 5 and len(self.contigs[0].seq) > self.levels[4].chunk_size and multipart_file:
        #     self.enable_fat_headers()  # first contig is huge and there's more contigs coming
        if len(self.contigs) > 10000:
            print("Over 10,000 scaffolds detected!  Titles for entries less than 10,000bp will not be drawn.")
            self.skip_small_titles = True
            self.sort_contigs = True  # Important! Skipping isn't valid unless they're sorted
        if self.sort_contigs:
            print("Scaffolds are being sorted by length.")
            self.contigs.sort(key=lambda fragment: -len(fragment.seq))  # Best to bring the largest contigs to the forefront

        for contig in self.contigs:  # Type: class DNASkittleUtils.Contigs.Contig
            length = len(contig.seq)
            title_length = len(contig.name) + 1  # for tracking where we are in the SEQUENCE file
            reset, title, tail = self.calc_padding(total_progress, length, multipart_file)

            contig.reset_padding = reset
            contig.title_padding = title
            contig.tail_padding = tail
            contig.nuc_title_start = seq_start
            contig.nuc_seq_start = seq_start + title_length

            total_progress += reset + title + tail + length  # pointer in image
            seq_start += title_length + length  # pointer in text
        return total_progress  # + reset + title + tail + length


    def read_contigs_and_calc_padding(self, input_file_path):
        try:
            self.contigs = read_contigs(input_file_path)
        except UnicodeDecodeError as e:
            print(e)
            print("Important: Non-standard characters detected.  Switching to 256 colormap for bytes")
            self.using_spectrum = True
            self.palette = viridis_palette()
            self.contigs = [Contig(input_file_path, open(input_file_path, 'rb').read())]
        return self.calc_all_padding()

    def prepare_image(self, image_length):
        width, height = self.max_dimensions(image_length)
        print("Image dimensions are", width, "x", height, "pixels")
        print("This will require approximately %s MB of RAM, or half that with --no_webpage" %
              "{:,}".format(width*height * 3 // 1048576 * 4))  # 3 channels, quadruple size for zoom tiles
        self.image = Image.new('RGB', (width, height), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.pixels = self.image.load()


    def calc_padding(self, total_progress, next_segment_length, multipart_file):
        min_gap = (20 + 6) * self.base_width  # 20px font height, + 6px vertical padding  * 100 nt per line
        if not multipart_file:
            return 0, 0, 0

        for i, current_level in enumerate(self.levels):
            if next_segment_length + min_gap < current_level.chunk_size:
                # give a full level of blank space just in case the previous
                title_padding = max(min_gap, self.levels[i - 1].chunk_size)
                if self.skip_small_titles and next_segment_length < small_title_bp:
                    # no small titles, but larger ones will still display,
                    title_padding = self.title_skip_padding  # normally 100 pixels per line
                    # this affects layout
                if not self.use_titles:
                    title_padding = 0  # don't leave space for a title, but still use tail and reset padding
                if title_padding > self.levels[3].chunk_size:  # Special case for full tile, don't need to go that big
                    title_padding = self.tile_label_size
                if next_segment_length + title_padding > current_level.chunk_size:
                    continue  # adding the title pushed the above comparison over the edge, step up one level
                space_remaining = current_level.chunk_size - total_progress % current_level.chunk_size
                # sequence comes right up to the edge.  There should always be >= 1 full gap
                reset_level = current_level  # bigger reset when close to filling chunk_size
                if next_segment_length + title_padding < space_remaining:
                    reset_level = self.levels[i - 1]
                # fill out the remainder so we can start at the beginning
                reset_padding = reset_level.chunk_size - total_progress % reset_level.chunk_size
                if total_progress == 0:  # nothing to reset from
                    reset_padding = 0
                total_padding = total_progress + title_padding + reset_padding + next_segment_length
                tail = self.levels[i - 1].chunk_size - total_padding % self.levels[i - 1].chunk_size - 1

                return reset_padding, title_padding, tail

        return 0, 0, 0


    def position_on_screen(self, index):
        """ Readable unoptimized version:
        Maps a nucleotide index to an x,y coordinate based on the rules set in self.levels"""
        xy = list(self.origin)  # column padding for various markup = self.levels[2].padding
        for i, level in enumerate(self.levels):
            if index < level.chunk_size:
                return int(xy[0]), int(xy[1])  # somehow a float snuck in here once
            part = i % 2
            coordinate_in_chunk = int(index / level.chunk_size) % level.modulo
            xy[part] += level.thickness * coordinate_in_chunk
        return int(xy[0]), int(xy[1])


    def draw_pixel(self, character, x, y):
        self.pixels[x, y] = self.palette[character]


    def draw_titles(self):
        total_progress = 0
        for contig in self.contigs:
            total_progress += contig.reset_padding  # is to move the cursor to the right line for a large title
            if contig.title_padding > self.title_skip_padding:  # there needs to be room to draw
                self.draw_title(total_progress, contig)
            total_progress += contig.title_padding + len(contig.seq) + contig.tail_padding


    def draw_title(self, total_progress, contig):
        upper_left = self.position_on_screen(total_progress)
        bottom_right = self.position_on_screen(total_progress + contig.title_padding - 2)
        width, height = bottom_right[0] - upper_left[0], bottom_right[1] - upper_left[1]

        font_size = 9  # font sizes: [9, 38, 380, 380 * 2]
        title_width = 18
        title_lines = 2

        # Title orientation and size
        vertical_label = contig.title_padding == self.levels[2].chunk_size
        if vertical_label:
            # column titles are vertically oriented
            width, height = height, width  # swap
            font_size = 38
            title_width = 50  # TODO: find width programatically
        if contig.title_padding >= self.levels[3].chunk_size:
            font_size = 380  # full row labels for chromosomes
            title_width = 50  # approximate width
        if contig.title_padding == self.tile_label_size:  # Tile dedicated to a Title (square shaped)
            # since this level is square, there's no point in rotating it
            font_size = 380 * 2  # doesn't really need to be 10x larger than the rows
            title_width = 50 // 2
            if self.use_fat_headers:
                # TODO add reset_padding from next contig, just in case there's unused space on this level
                tiles_spanned = int(math.ceil((len(contig.seq) + contig.tail_padding) / self.levels[4].chunk_size))
                title_width *= tiles_spanned  # Twice the size, but you have 3 tile columns to fill, also limited by 'width'
                title_lines = 1
                upper_left[1] -= self.levels[3].thickness  # above the start of the coordinate grid
                height = self.levels[3].thickness
                width = self.levels[4].thickness * tiles_spanned  # spans 3 full Tiles, or one full Page width

        contig_name = contig.name
        self.write_title(contig_name, width, height, font_size, title_lines, title_width, upper_left, vertical_label)


    def write_title(self, contig_name, width, height, font_size, title_lines, title_width, upper_left, vertical_label):
        upper_left = list(upper_left)  # to make it mutable
        if font_size in self.fonts:
            font = self.fonts[font_size]
        else:
            font = ImageFont.truetype(font_name, font_size)
        multi_line_title = pretty_contig_name(contig_name, title_width, title_lines)
        txt = Image.new('RGBA', (width, height))
        bottom_justified = height - multi_line_height(font, multi_line_title, txt)
        ImageDraw.Draw(txt).multiline_text((0, max(0, bottom_justified)), multi_line_title, font=font,
                                           fill=(0, 0, 0, 255))
        if vertical_label:
            txt = txt.rotate(90, expand=True)
            upper_left[0] += 8  # adjusts baseline for more polish
        self.image.paste(txt, (upper_left[0], upper_left[1]), txt)


    def output_image(self, output_folder, output_file_name):
        del self.pixels
        del self.draw
        self.final_output_location = os.path.join(output_folder, output_file_name + ".png")
        print("-- Writing:", self.final_output_location, "--")
        self.image.save(self.final_output_location, 'PNG')
        # del self.image


    def max_dimensions(self, image_length):
        """ Uses Tile Layout to find the largest chunk size in each dimension (XY) that the
        image_length will reach
        :param image_length: includes sequence length and padding from self.read_contigs_and_calc_padding()
        :return: width and height needed
        """
        width_height = [0, 0]
        for i, level in enumerate(self.levels):
            part = i % 2
            # how many of these will you need up to a full modulo worth
            coordinate_in_chunk = min(int(math.ceil(image_length / float(level.chunk_size))), level.modulo)
            if coordinate_in_chunk > 1:
                # not cumulative, just take the max size for either x or y
                width_height[part] = max(width_height[part], level.thickness * coordinate_in_chunk)
        width_height = [sum(x) for x in zip(width_height, self.origin)]  # , [self.levels[2].padding] * 2
        width_height[0] += self.levels[2].padding   # add column padding to both sides
        width_height[1] += self.levels[2].padding   # column padding used as a proxy for vertical padding
        width_height[0] += self.origin[0]  # add in origin offset
        width_height[1] += self.origin[1]
        return int(width_height[0]), int(width_height[1])


    def generate_html(self, output_folder, output_file_name):
        try:
            import DDV
            module_path = os.path.dirname(DDV.__file__)
            html_template = os.path.join(module_path, 'html_template')
            copytree(html_template, output_folder)  # copies the whole template directory
            html_path = os.path.join(output_folder, 'index.html')
            html_content = {"title": output_file_name.replace('_', ' '),
                            "originalImageWidth": str(self.image.width if self.image else 1),
                            "originalImageHeight": str(self.image.height if self.image else 1),
                            "image_origin": str(self.origin),
                            "ColumnPadding": str(self.levels[2].padding),
                            "columnWidthInNucleotides": str(self.levels[1].chunk_size),
                            "layoutSelector": '1',
                            "layout_levels": self.levels_json(),
                            "ContigSpacingJSON": self.contig_json(),
                            "multipart_file": str(len(self.contigs) > 1).lower(),
                            # "use_fat_headers": str(self.use_fat_headers).lower(),  # use image_origin and layout_levels
                            "includeDensity": 'false',
                            "ipTotal": str(self.image_length),
                            "direct_data_file_length": str(self.image_length),  # TODO: this isn't right because includes padding
                            "sbegin": '1',
                            "send": str(self.image_length),
                            "date": datetime.now().strftime("%Y-%m-%d")}
            html_content['legend'] = """    <strong>Legend:</strong>
                <img class='legend-icon' src='img/LEGEND-A.png'/>
                <img class='legend-icon' src='img/LEGEND-T.png'/>
                <img class='legend-icon' src='img/LEGEND-G.png'/>
                <img class='legend-icon' src='img/LEGEND-C.png'/>
                <img class='legend-icon' src='img/LEGEND-N.png'/>
                <img class='legend-icon' src='img/LEGEND-bg.png'/>
                <span class='color-explanation'>Color blind safe colors.  G/C rich regions are blue/green.
                    A/T rich areas are reddish.  Poly-purines are more yellow (orange/green).
                    Poly-pyrimidines are more purple (blue/red).
                    Diffuse natural colors were chosen to be less harsh on the eyes.</span>
            """
            if not self.low_contrast:
                html_content['legend'] = """    <strong>Legend:</strong>
                                <img class='legend-icon' src='img/LEGEND-A-contrast.png'/>
                                <img class='legend-icon' src='img/LEGEND-T-contrast.png'/>
                                <img class='legend-icon' src='img/LEGEND-G-contrast.png'/>
                                <img class='legend-icon' src='img/LEGEND-C-contrast.png'/>
                                <img class='legend-icon' src='img/LEGEND-N.png'/>
                                <img class='legend-icon' src='img/LEGEND-bg.png'/>
                                <span class='color-explanation'>G/C rich regions are blue/green.
                                    A/T rich areas are reddish.</span>
                            """
            if self.using_spectrum:
                html_content['legend'] = """    <strong>Legend:</strong>
                <span class='color-explanation'>Each pixel is 1 byte with a range of 0 - 255. 
                0 = dark purple. 125 = green, 255 = yellow. Developed as 
                Matplotlib's default color palette.  It is 
                perceptually uniform and color blind safe.</span>"""
            html_content.update(self.additional_html_content(html_content))
            with open(os.path.join(html_template, 'index.html'), 'r') as template:
                template_content = template.read()
                for key, value in html_content.items():
                    template_content = template_content.replace('{{' + key + '}}', value)
                with open(html_path, 'w') as out:
                    out.write(template_content)

        except Exception as e:
            print('Error while generating HTML:', '\n')
            traceback.print_exc()


    def contig_struct(self):
        json = []
        xy_seq_start = 0
        for index, contig in enumerate(self.contigs):
            if index > 1000:
                break  # I don't want to use a slice operator on the for loop because that will copy it
            xy_seq_start += contig.reset_padding + contig.title_padding
            xy_seq_end = xy_seq_start + len(contig.seq)
            json.append({"name": contig.name.replace("'", ""), "xy_seq_start": xy_seq_start, "xy_seq_end": xy_seq_end,
                         "title_padding": contig.title_padding, "tail_padding": contig.tail_padding,
                         "xy_title_start": xy_seq_start - contig.title_padding,
                         "nuc_title_start": contig.nuc_title_start, "nuc_seq_start": contig.nuc_seq_start})
            xy_seq_start += len(contig.seq) + contig.tail_padding
        return json

    def contig_json(self):
        json = self.contig_struct()
        return "[" + ',\n'.join([str(x) for x in json]) + "]"


    def levels_json(self):
        json = []
        for level in self.levels:
            json.append({"modulo": level.modulo, "chunk_size": level.chunk_size,
                         "padding": level.padding, "thickness": level.thickness})
        return str(json)

    def get_packed_coordinates(self):
        """An attempted speed up for draw_nucleotides() that was the same speed.  In draw_nucleotides() the
        extra code was:
            coordinates = self.get_packed_coordinates()  # precomputed sets of 100,000
            seq_consumed = 0
            columns_batched = 0
            for column in range(0, seq_length, 100000):
                if seq_length - column > 100000:
                    columns_batched += 1
                    x, y = self.position_on_screen(total_progress)  # only one call per column
                    for cx, cy, offset in coordinates:
                        self.draw_pixel(contig.seq[column + offset], x + cx, y + cy)
                    total_progress += 100000
                    seq_consumed += 100000
                else:
                    pass  # loop will exit and remaining seq will be handled individually

        This method is an optimization that computes all offsets for a column once so they can be reused.
        The output looks like this:  (x, y, sequence offset)
        [(0, 0, 0), (1, 0, 1), (2, 0, 2), (3, 0, 3), ... (0, 1, 10), (1, 1, 11), (2, 1, 12), (3, 1, 13),"""
        line = range(self.levels[0].modulo)
        column_height = self.levels[1].modulo
        coords = []
        for y in range(column_height):
            coords.extend([(x, y, y * self.levels[0].modulo + x) for x in line])
        return coords


    def additional_html_content(self, html_content):
        return {}  # override in children

