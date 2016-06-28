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
import math
import textwrap

from PIL import Image, ImageDraw, ImageFont
from datetime import datetime
from collections import defaultdict
from deepzoom import deepzoom

# Original DDV Colors
palette = defaultdict(lambda: (0, 0, 0))
palette['A'] = (255, 0, 0)
palette['a'] = (255, 0, 0)  # A
palette['G'] = (0, 255, 0)
palette['g'] = (0, 255, 0)  # G
palette['T'] = (250, 240, 114)
palette['t'] = (250, 240, 114)  # T
palette['C'] = (0, 0, 255)
palette['c'] = (0, 0, 255)  # C
palette['N'] = (30, 30, 30)
palette['n'] = (30, 30, 30)  # N


class LayoutLevel:
    def __init__(self, name, modulo, chunk_size=None, padding=0, thickness=1, levels=None):
        self.modulo = modulo
        if chunk_size is not None:
            self.chunk_size = chunk_size
            self.padding = padding
            self.thickness = thickness
        else:
            child = levels[-1]
            self.chunk_size = child.modulo * child.chunk_size
            self.padding = 6 * int(3 ** (len(levels) - 2))  # third level (count=2) should be 6, then 18
            last_parallel = levels[-2]
            self.thickness = last_parallel.modulo * last_parallel.thickness + self.padding


class Contig:
    def __init__(self, name, seq, reset_padding, title_padding, tail_padding):
        self.name = name
        self.seq = seq
        self.reset_padding = reset_padding
        self.title_padding = title_padding
        self.tail_padding = tail_padding


class DDVTileLayout:
    def __init__(self):
        # use_fat_headers: For large chromosomes in multipart files, do you change the layout to allow for titles that
        # are outside of the nucleotide coordinate grid?
        self.use_fat_headers = False  # Can only be changed in code.
        self.image = None
        self.draw = None
        self.pixels = None
        self.contigs = []
        self.image_length = 0
        # noinspection PyListCreation
        self.levels = [
            LayoutLevel("XInColumn", 100, 1),  # [0]
            LayoutLevel("LineInColumn", 1000, 100)  # [1]
        ]
        self.levels.append(LayoutLevel("ColumnInRow", 100, levels=self.levels))  # [2]
        self.levels.append(LayoutLevel("RowInTile", 10, levels=self.levels))  # [3]
        self.levels.append(LayoutLevel("XInTile", 3, levels=self.levels))  # [4]
        self.levels.append(LayoutLevel("YInTile", 4, levels=self.levels))  # [5]
        if self.use_fat_headers:
            self.levels[5].padding += self.levels[3].thickness  # one full row for a chromosome title
            self.levels[5].thickness += self.levels[3].thickness
        self.levels.append(LayoutLevel("TileColumn", 9, levels=self.levels))  # [6]
        self.levels.append(LayoutLevel("TileRow", 999, levels=self.levels))  # [7]

        self.tile_label_size = self.levels[3].chunk_size * 2
        if self.use_fat_headers:
            self.tile_label_size = 0  # Fat_headers are not part of the coordinate space


    def process_file(self, input_file_name, output_folder, output_file_name):
        start_time = datetime.now()
        self.image_length = self.read_contigs(input_file_name)
        print("Read contigs :", datetime.now() - start_time)
        self.prepare_image(self.image_length)
        print("Initialized Image:", datetime.now() - start_time, "\n")
        try:  # These try catch statements ensure we get at least some output.  These jobs can take hours
            self.draw_nucleotides()
            print("\nDrew Nucleotides:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing nucleotides:', '\n', str(e))
        try:
            if len(self.contigs) > 1:
                self.draw_titles()
                print("Drew Titles:", datetime.now() - start_time)
        except Exception as e:
            print('Encountered exception while drawing titles:', '\n', str(e))
        try:
            self.generate_html(input_file_name, output_folder, output_file_name)
        except: pass
        self.output_image(output_folder, output_file_name)
        print("Output Image in:", datetime.now() - start_time)

    def draw_nucleotides(self):
        total_progress = 0
        # if self.image_length > 300000000:
        #     positioner = self.position_on_screen_big  # big images
        # else:
        #     positioner = self.position_on_screen  # small images
        # Layout contigs one at a time
        for contig in self.contigs:
            total_progress += contig.reset_padding + contig.title_padding
            # worker.ReportProgress((int) (nucleotidesProcessed += contig.len(seq))) # doesn't include padding
            for cx in range(0, len(contig.seq), 100):
                x, y = self.position_on_screen(total_progress)
                remaining = min(100, len(contig.seq) - cx)
                total_progress += remaining
                for i in range(remaining):
                    self.draw_pixel(contig.seq[cx + i], x + i, y)
            if self.image_length > 10000000:
                print('\r', str(total_progress / self.image_length * 100)[:6], '% done:', contig.name,
                      end="")  # pseudo progress bar
            total_progress += contig.tail_padding  # add trailing white space after the contig sequence body
        print()

    def read_contigs(self, input_file_name):
        self.contigs = []
        total_progress = 0
        current_name = ""
        seq_collection = []

        # Pre-read generates an array of contigs with labels and sequences
        with open(input_file_name, 'r') as streamFASTAFile:
            for read in streamFASTAFile.read().splitlines():
                if read == "":
                    continue
                if read[0] == ">":
                    # If we have sequence gathered and we run into a second (or more) block
                    if len(seq_collection) > 0:
                        sequence = "".join(seq_collection)
                        seq_collection = []  # clear
                        reset, title, tail = self.calc_padding(total_progress, len(sequence), True)
                        self.contigs.append(Contig(current_name, sequence, reset, title, tail))
                        total_progress += reset + title + tail + len(sequence)
                    current_name = read[1:]  # remove >
                else:
                    # collects the sequence to be stored in the contig, constant time performance don't concat strings!
                    seq_collection.append(read)

        # add the last contig to the list
        sequence = "".join(seq_collection)
        reset, title, tail = self.calc_padding(total_progress, len(sequence), len(self.contigs) > 0)
        self.contigs.append(Contig(current_name, sequence, reset, title, tail))
        return total_progress + reset + title + tail + len(sequence)


    def prepare_image(self, image_length):
        width, height = self.max_dimensions(image_length)
        self.image = Image.new('RGB', (width, height), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.pixels = self.image.load()


    def calc_padding(self, total_progress, next_segment_length, multipart_file):
        min_gap = (20 + 6) * 100  # 20px font height, + 6px vertical padding  * 100 nt per line
        if not multipart_file:
            return 0, 0, 0
        
        for i, current_level in enumerate(self.levels):
            if next_segment_length + min_gap < current_level.chunk_size:
                # give a full level of blank space just in case the previous
                title_padding = max(min_gap, self.levels[i - 1].chunk_size)
                if title_padding > self.levels[3].chunk_size:  # Special case for full tile, don't need to go that big
                    title_padding = self.tile_label_size
                if next_segment_length + title_padding > current_level.chunk_size:
                    continue  # adding the title pushed the above comparison over the edge, step up one level
                space_remaining = current_level.chunk_size - total_progress % current_level.chunk_size
                # sequence comes right up to the edge.  There should always be >= 1 full gap
                reset_level = current_level  # bigger reset when close to filling chunk_size
                if next_segment_length + title_padding < space_remaining:
                    reset_level = self.levels[i - 1]
                    # reset_level = self.levels[i - 1]
                reset_padding = 0
                if total_progress != 0:  # fill out the remainder so we can start at the beginning
                    reset_padding = reset_level.chunk_size - total_progress % reset_level.chunk_size
                total_padding = total_progress + title_padding + reset_padding + next_segment_length
                tail = self.levels[i - 1].chunk_size - total_padding % self.levels[i - 1].chunk_size - 1

                return reset_padding, title_padding, tail

        return 0, 0, 0


    def position_on_screen(self, index):
        """ Readable unoptimized version:
        Maps a nucleotide index to an x,y coordinate based on the rules set in self.levels"""
        xy = [0, 0]
        if self.use_fat_headers:
            xy[1] = self.levels[5].padding  # padding comes before, not after
        for i, level in enumerate(self.levels):
            if index < level.chunk_size:
                return xy
            part = i % 2
            coordinate_in_chunk = int(index / level.chunk_size) % level.modulo
            xy[part] += level.thickness * coordinate_in_chunk
        return xy

    @staticmethod
    def position_on_screen_small(index):
        # Less readable
        # x = self.levels[0].thickness * (int(index / self.levels[0].chunk_size) % self.levels[0].modulo)
        # x+= self.levels[2].thickness * (int(index / self.levels[2].chunk_size) % self.levels[2].modulo)
        # y = self.levels[1].thickness * (int(index / self.levels[1].chunk_size) % self.levels[1].modulo)
        # y+= self.levels[3].thickness * (int(index / self.levels[3].chunk_size) % self.levels[3].modulo)

        x = index % 100 + 106 * ((index // 100000) % 100) + 10654 * (index // 100000000)  # % 3)
        y = (index // 100) % 1000 + 1018 * ((index // 10000000) % 10)  # + 10342 * ((index // 300000000) % 4)
        return x, y

    @staticmethod
    def position_on_screen_big(index):
        # 10654 * 3 + 486 padding = 32448
        x = index % 100 + 106 * ((index // 100000) % 100) + 10654 * ((index // 100000000) % 3) + \
            32448 * (index // 1200000000)  # % 9 #this will continue tile columns indefinitely (8 needed 4 human genome)
        y = (index // 100) % 1000 + 1018 * ((index // 10000000) % 10) + 10342 * ((index // 300000000) % 4)
        # + 42826 * (index // 10800000000)  # 10342 * 4 + 1458 padding = 42826
        return [x, y]

    def draw_pixel(self, character, x, y):
        self.pixels[x, y] = palette[character]


    def draw_titles(self):
        total_progress = 0
        for contig in self.contigs:
            total_progress += contig.reset_padding  # is to move the cursor to the right line for a large title
            self.draw_title(total_progress, contig)
            total_progress += contig.title_padding + len(contig.seq) + contig.tail_padding

    def draw_title(self, total_progress, contig):
        upper_left = self.position_on_screen(total_progress)
        bottom_right = self.position_on_screen(total_progress + contig.title_padding - 2)
        width, height = bottom_right[0] - upper_left[0], bottom_right[1] - upper_left[1]

        font_size = 9
        title_width = 18
        title_lines = 2

        # Title orientation and size
        if contig.title_padding == self.levels[2].chunk_size:
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
                tiles_spanned = math.ceil((len(contig.seq) + contig.tail_padding) / self.levels[4].chunk_size)
                title_width *= tiles_spanned  # Twice the size, but you have 3 tile columns to fill, also limited by 'width'
                title_lines = 1
                upper_left[1] -= self.levels[3].thickness  # above the start of the coordinate grid
                height = self.levels[3].thickness
                width = self.levels[4].thickness * tiles_spanned  # spans 3 full Tiles, or one full Page width

        font = ImageFont.truetype("tahoma.ttf", font_size)
        txt = Image.new('RGBA', (width, height))
        multi_line_title = pretty_contig_name(contig, title_width, title_lines)

        bottom_justified = height - multi_line_height(font, multi_line_title, txt)
        ImageDraw.Draw(txt).multiline_text((0, max(0, bottom_justified)), multi_line_title, font=font, fill=(0, 0, 0, 255))
        if contig.title_padding == self.levels[2].chunk_size:
            txt = txt.rotate(90, expand=True)
            upper_left[0] += 8  # adjusts baseline for more polish

        self.image.paste(txt, (upper_left[0], upper_left[1]), txt)

    def output_image(self, output_folder, output_file_name):
        del self.pixels
        del self.draw
        self.image.save(os.path.join(output_folder, output_file_name), 'PNG')
        del self.image


    def max_dimensions(self, image_length):
        """ Uses Tile Layout to find the largest chunk size in each dimension (XY) that the
        image_length will reach
        :param image_length: includes sequence length and padding from self.read_contigs()
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
        if self.use_fat_headers:  # extra margin at the top of the image for a title
            width_height[1] += self.levels[5].padding
        return width_height

    def generate_html(self, input_file_name, output_folder, output_file_name):
        os.makedirs(output_folder, exist_ok=True)
        html_path = os.path.join(output_folder, 'embed.html')
        html_content = """
        <!DOCTYPE html>
        <html lang='en'>
        <head>
        <meta http-equiv='Content-Type' content='text/html; charset=UTF-8' />
        <title>DNA Data Visualization : """ + input_file_name[:input_file_name.rfind('.')] + """</title>
        <script src='../../openseadragon.js' type='text/javascript'></script>
        <script type='text/javascript' src='../../jquery-1.7.min.js'></script>
        <script src='../../openseadragon-scalebar.js' type='text/javascript'></script>

        <script type='text/javascript'>
	        var originalImageWidth= """ + str(self.image.width) + """;
            var originalImageHeight= """ + str(self.image.height) + """;
            var ColumnPadding = """ + str(self.levels[2].padding) + """;
            var columnWidthInNucleotides = """ + str(self.levels[1].chunk_size) + """;
            var layoutSelector = 1;
            var layout_levels = """ + self.levels_json() + """;
            var ContigSpacingJSON = """ + self.contig_json() + """;
            var multipart_file = """ + str(len(self.contigs) > 1).lower() + """;
            var use_fat_headers = """ + str(self.use_fat_headers).lower() + """;
            var includeDensity = false;

            var usa='refseq_fetch:""" + input_file_name + """';
            var ipTotal = """ + str(self.image_length) + """;
            var direct_data_file='sequence.fasta';
            var direct_data_file_length=1;
            var sbegin='1';
            var send=""" + str(self.image_length) + """;
            var output_dir = '../../';
            </script>
        <script src='../../nucleotideNumber.js' type='text/javascript'></script>
        <script src='../../nucleicDensity.js' type='text/javascript'></script>
        <link rel='stylesheet' type='text/css' href='../../seadragon.css' />
            <!-- BIOJS css -->
            <script language='JavaScript' type='text/javascript' src='../../Biojs.js'></script>
            <!-- component code -->
            <script language='JavaScript' type='text/javascript' src='../../Biojs.Sequence.js'></script>
        </head>

        <body>
        <h2 class='mainTitle'>Data Visualization - DNA</h2>
        <span style='float:left;'>Menu:&nbsp;</span>
            <ul class='selectChromosome'>
            <li><a href='../'>Select Visualization</a></li>
             </ul>
        <h2 class='mainTitle'><strong>""" + input_file_name + """</strong>
         </h2>

        <div id='container' class='chromosome-container' data-chr-source=''>
        </div>

        <p class='legendHeading'><strong>Legend:</strong><br /></p>
        <div style='margin-left:50px;'>
            <img src='../../LEGEND-A.png' />
            <img src='../../LEGEND-T.png' />
            <img src='../../LEGEND-G.png' />
            <img src='../../LEGEND-C.png' />
            <img src='../../LEGEND-N.png' />
            <img src='../../LEGEND-bg.png' />
        </div>

        <script type='text/javascript'>
            outputTable();
            if (includeDensity) { outputDensityUI();}
        </script>

        <div class='legend-details'>
        <h3>Data Source:</h3>
        <a href='sequence.fasta'>Download FASTA file</a>
        <br />
        Custom/local sequence (DDV seq ID): """ + input_file_name + """<br />

        <h3>Notes</h3>
        This DNA data visualization interface was generated with <a href='https://github.com/photomedia/DDV'>DDV</a>
        <br />Date Visualization Created:""" + datetime.now().strftime("%Y-%m-%d") + """
        <script type='text/javascript'>
                    otherCredits();
        </script>
        </div>
        </body>

        </html>
        """
        with open(html_path, 'w') as out:
            out.write(html_content)

    def contig_json(self):
        json = []
        startingIndex = 0  # camel case is left over from C# for javascript compatibility
        for contig in self.contigs:
            startingIndex += contig.title_padding
            endIndex = startingIndex + len(contig.seq)
            json.append({"name": contig.name.replace("'", ""), "startingIndex": startingIndex, "endIndex": endIndex,
                         "title_padding": contig.title_padding, "tail_padding": contig.tail_padding})
        return str(json)

    def levels_json(self):
        json = []
        for level in self.levels:
            json.append({"modulo": level.modulo, "chunk_size": level.chunk_size,
                         "padding": level.padding, "thickness": level.thickness})
        return str(json)


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


def create_deepzoom_stack(input_image, output_dzi):
    dz_params = {'tile_size': 256,
                               'tile_overlap': 1,
                               'tile_format': "jpg",
                               'image_quality': 0.85,
                               'resize_filter': "antialias"}

    creator = deepzoom.ImageCreator(tile_size=dz_params['tile_size'],
                                    tile_overlap=dz_params['tile_overlap'],
                                    tile_format=dz_params['tile_format'],
                                    image_quality=dz_params['image_quality'],
                                    resize_filter=dz_params['resize_filter'])

    creator.create(input_image, output_dzi)


if __name__ == "__main__":
    from ParallelGenomeLayout import ParallelLayout

    # input_file_name, output_file_name = "sequence.fa", "output.png"

    # layout.process_file("Animalia_Mammalia_Homo_Sapiens_GRCH38_chr20.fa", ".", "ch20-2.png")
    # layout.process_file("Animalia_Mammalia_Homo_Sapiens_GRCH38_nonchromosomal.fa", ".", "non-chromosomal.png")
    # layout.process_file("Animalia_Mammalia_Homo_Sapiens_GRCH38_chr1.fa", ".", "chr1 Human.png")
    # layout.process_file("Human selenoproteins.fa", ".", "selenoproteins.png")
    # layout.process_file("multi_part_layout.fa", ".", "multi_part_layout.png")
    # layout.process_file("CYUI01000001-CYUI01015997.fasta", ".", "susie3.png")
    folder = "."
    image = sys.argv[1][:sys.argv[1].rfind(".")] + ".png"
    n_arguments = len(sys.argv)
    if n_arguments >= 4:
        folder, image = sys.argv[2], sys.argv[3]
    if n_arguments == 3:
        image = sys.argv[2]

    if n_arguments == 2:
        output_file = os.path.basename(sys.argv[1].replace('.png', '.dzi'))
        output_dir = os.path.dirname(sys.argv[1])

        create_deepzoom_stack(sys.argv[1], os.path.join(output_dir, output_file))
        sys.exit(0)
    elif n_arguments > 4:  # Parallel genome column layout
        layout = ParallelLayout(n_arguments - 3)
        layout.process_file(sys.argv[1], folder, image, sys.argv[4:])
    else:
        layout = DDVTileLayout()
        layout.process_file(sys.argv[1], folder, image)

    create_deepzoom_stack(os.path.join(folder, image), os.path.join(folder, str(image).replace('.png', '.dzi')))
    sys.exit(0)
