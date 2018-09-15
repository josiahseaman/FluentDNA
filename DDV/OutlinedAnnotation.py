from itertools import chain

from PIL import Image, ImageFont, ImageDraw

from DDV.Annotations import GFF, extract_gene_name, find_universal_prefix
from DDV.DDVUtils import multi_line_height
from DDV.Span import Span
from DDV.TileLayout import TileLayout, hex_to_rgb
from DDV.DDVUtils import linspace
from collections import namedtuple
Point = namedtuple('Point', ['x', 'y'])


def blend_pixel(markup_canvas, pt, c, overwrite=False):
    if overwrite or markup_canvas[pt[0], pt[1]][3] == 0:  # nothing drawn
        markup_canvas[pt[0], pt[1]] = c
    else:
        remaining_light = 1.0 - (markup_canvas[pt[0], pt[1]][3] / 256)
        combined_alpha = 256 - int(remaining_light * (256 - c[3]) )
        markup_canvas[pt[0], pt[1]] = (c[0], c[1], c[2], combined_alpha)

def annotation_points(entry, renderer, progress_offset):
    annotation_points = []
    for i in range(entry.start, entry.end):
        # important to include title and reset padding in coordinate frame
        progress = i + progress_offset
        annotation_points.append(renderer.position_on_screen(progress))
    return annotation_points


class OutlinedAnnotation(TileLayout):
    def __init__(self, gff_file, query=None, **kwargs):
        super(OutlinedAnnotation, self).__init__(**kwargs)
        self.annotation = GFF(gff_file).annotations if gff_file is not None else None
        self.query_annotation = GFF(query).annotations if query is not None else None
        self.pil_mode = 'RGBA'  # Alpha channel necessary for outline blending
        self.font_name = "ariblk.ttf"  # TODO: compatibility testing with Mac
        self.border_width = 30

    def process_file(self, input_file_path, output_folder, output_file_name,
                     no_webpage=False, extract_contigs=None):
        if self.annotation is not None:
            with open(input_file_path, 'r') as fasta:
                assert fasta.readline().startswith('>'), "Fasta file must start with a header '>name'"
        super(OutlinedAnnotation, self).process_file(input_file_path, output_folder, output_file_name,
                                                     no_webpage, extract_contigs)
        # nothing extra

    def draw_extras(self):
        """Drawing Annotations labels and shadow outlines"""
        markup_image = Image.new('RGBA', (self.image.width, self.image.height), (0,0,0,0))
        markup_canvas = markup_image.load()
        annotated_regions, query_regions = [], []

        if self.annotation is not None:
            annotated_regions = self.draw_annotation_outlines(self.annotation, markup_canvas, (65, 42, 80))
        if self.query_annotation is not None:
            query_regions = self.draw_annotation_outlines(self.query_annotation, markup_canvas, (211, 229, 199))
        # important to combine set of names so not too much prefix gets chopped off
        annotated_regions = list(chain(annotated_regions, query_regions))
        universal_prefix = find_universal_prefix(annotated_regions)
        print("Removing Universal Prefix from annotations:", universal_prefix)
        self.draw_annotation_labels(markup_image, annotated_regions, universal_prefix)
        self.image = Image.alpha_composite(self.image, markup_image)

    def draw_annotation_outlines(self, annotations, markup_canvas, shadow_color):
        regions = self.find_annotated_regions(annotations)
        self.draw_exons(markup_canvas, regions)
        annotation_point_union = self.draw_big_shadow_outline(markup_canvas, regions, shadow_color)
        self.draw_secondary_shadows(annotation_point_union, markup_canvas, regions, shadow_color)
        return regions

    def draw_big_shadow_outline(self, markup_canvas, regions, shadow):
        print("Drawing annotation outlines")
        annotation_point_union = set()
        for region in regions:
            annotation_point_union.update(region.points)
            region.outline_points = outlines(region.points,  # small outline
                                             self.border_width // 4, self.image.width, self.image.height)
        # desaturated purple drop shadow, decreasing opacity
        opacities = linspace(197, 10, self.border_width)
        outline_colors = [(shadow[0], shadow[1], shadow[2], int(opacity)) for opacity in opacities]
        big_shadow = outlines(annotation_point_union,
                              self.border_width, self.image.width, self.image.height)
        self.draw_shadow(big_shadow, markup_canvas, outline_colors)
        return annotation_point_union

    def draw_exons(self, markup_canvas, regions):
        print("Drawing exons")
        exon_color = (255, 255, 255, 107)  # white highlighter.  This is less disruptive overall
        for region in regions:
            for point in region.exon_region_points():  # highlight exons
                blend_pixel(markup_canvas, point, exon_color)

    def draw_secondary_shadows(self, annotation_point_union, markup_canvas, regions, shadow):
        """Find subset of genes who are completely overshadowed"""
        print("Drawing secondary shadows")
        opacities = linspace(170, 40, self.border_width // 4)
        outline_colors = [(shadow[0], shadow[1], shadow[2], int(opacity)) for opacity in opacities]
        for region in regions:
            # if annotation_point_union.issuperset(region.outline_points):
            lost_edge = region.outline_points[0].intersection(annotation_point_union)
            if lost_edge:
                lost_shadows = [layer.intersection(annotation_point_union) for layer in region.outline_points]
                self.draw_shadow(lost_shadows, markup_canvas, outline_colors)

    def draw_shadow(self, shadow, markup_canvas, outline_colors, flat_color=False):
        for radius, layer in enumerate(shadow):
            darkness = radius
            # self.border_width - len(region.outline_points) + radius  # softer line for small features
            c = outline_colors[darkness]
            for pt in layer:
                blend_pixel(markup_canvas, pt, c)

    def find_annotated_regions(self, annotations):
        """ :type annotations: dict(GFF.Annotation) """
        print("Collecting points in annotated regions")
        positions = self.contig_struct()
        regions = []

        for sc_index, coordinate_frame in enumerate(positions):  # Exact match required (case sensitive)
            genes_seen = set()  # keeping track so I don't double down on gene/mRNA/transcript
            scaff_name = coordinate_frame["name"].split()
            if not scaff_name:
                raise ValueError('Annotation cannot proceed without a contig name in the FASTA file.  \n'
                                 'Please add a line at the beginning of your fasta file with a name '
                                 'that exactly matches the first column in your annotation. For example: '
                                 '>chrMt')
            scaff_name = scaff_name[0]
            if scaff_name in annotations.keys():
                for entry in annotations[scaff_name]:
                    # redundancy checks for file with both mRNA and gene
                    if entry.feature == 'gene':
                        regions.append(AnnotatedRegion(entry, self, coordinate_frame["xy_seq_start"]))
                        genes_seen.add(extract_gene_name(entry).replace('g','t'))
                    if entry.feature == 'mRNA' and extract_gene_name(entry) not in genes_seen:
                        regions.append(AnnotatedRegion(entry, self, coordinate_frame["xy_seq_start"]))
                    if entry.feature == 'CDS' or entry.feature == 'exon':
                        # hopefully mRNA comes first in the file
                        # if extract_gene_name(regions[-1]) == entry.attributes['Parent']:
                        # the if was commented out because CDS [Parent] to mRNA, not gene names
                            regions[-1].add_cds_region(entry)

        return regions


    def draw_annotation_labels(self, markup_image, annotated_regions, universal_prefix=''):
        """ :type annotated_regions: list(AnnotatedRegion) """
        print("Drawing annotation labels")
        self.fonts = {9: ImageFont.load_default()}  # clear font cache, this may be a different font
        for region in annotated_regions:
            try:
                pts = [pt for pt in region.points]
                left, right = min(pts, key=lambda p: p[0])[0], max(pts, key=lambda p: p[0])[0]
                top, bottom = min(pts, key=lambda p: p[1])[1], max(pts, key=lambda p: p[1])[1]

                width, height, left, right, top = self.handle_multi_column_annotations(region, left, right,
                                                                                       top, bottom)
                vertical_label = height > width
                upper_left = [left, top]

                # Title orientation and size
                if vertical_label:
                    width, height = height, width  # swap

                font_size_by_width  = max(9, int((min(3000, width) * 0.09)))  # found eq with two reference points
                font_size_by_height = max(9, int((min(3000, height * 18) * 0.09)))
                if height <= 244: # 398580bp = 1900 width, 243 height
                    font_size_by_height = min(font_size_by_height, int(1900 * .09))  # 171 max font size in one fiber
                font_size = min(font_size_by_width, font_size_by_height)
                if height < 11:
                    height = 11  # don't make the area so small it clips the text
                    upper_left[1] -= 2

                self.write_label(extract_gene_name(region, universal_prefix), width, height, font_size, 18,
                                 upper_left, vertical_label, region.strand, markup_image)
            except BaseException as e:
                print('Error while drawing label %s' % extract_gene_name(region), e)

    def handle_multi_column_annotations(self, region, left, right, top, bottom):
        multi_column = abs(right - left) > self.base_width
        if multi_column:  # pick the biggest column to contain the label, ignore others
            median_point = len(region.points) // 2 + region.xy_seq_start + min(region.start, region.end)
            s = median_point // self.base_width * self.base_width  # beginning of the line holding median
            left = self.position_on_screen(s)[0]  # x coordinate of beginning of line
            right = self.position_on_screen(s + self.base_width - 1)[0]  # end of one line
            filtered = [pt for pt in region.points if right > pt[0] > left]  # off by ones here don't matter
            top, bottom = min(filtered, key=lambda p: p[1])[1], max(filtered, key=lambda p: p[1])[1]
            height = len(filtered) // self.base_width
        else:
            height = len(region.points) // self.base_width
        width = self.base_width
        return width, height, left, right, top

    def write_label(self, contig_name, width, height, font_size, title_width, upper_left, vertical_label,
                    strand, canvas, horizontal_centering=False, center_vertical=False, chop_text=True):
        """write_label() made to nicely draw single line gene labels from annotation
        :param horizontal_centering:
        """
        font = self.get_font(self.font_name, font_size)
        upper_left = list(upper_left)  # to make it mutable
        shortened = contig_name[-title_width:]  # max length 18.  Last characters are most unique
        txt = Image.new('RGBA', (width, height))
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
        text_color = (0, 0, 0, 255) if font_size < 14 else (50, 50, 50, 235)
        if font_size > 30:
            text_color = (100, 100, 100, 200)
        txt_canvas.multiline_text((0, max(0, vertically_centered)), shortened, font=font,
                                           fill=text_color)
        if vertical_label:
            rotation_direction = 90 if strand == '-' else -90
            txt = txt.rotate(rotation_direction, expand=True)
            upper_left[1] += -4 if strand == '-' else 4
        if horizontal_centering:
            margin = width - text_width
            upper_left[0] += margin // 2
        canvas.paste(txt, (upper_left[0], upper_left[1]), txt)


def getNeighbors(x, y):
    return (x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)

def allNeighbors(x, y):
    return set(getNeighbors(x, y)).union({(x + 1, y + 1), (x - 1, y - 1),
                                   (x - 1, y + 1), (x + 1, y - 1)})

def outlines(annotation_points, radius, width, height):
    workingSet = set(annotation_points)
    nextEdge = set(annotation_points)
    layers = []
    for iterationStep in range(radius, 0,  -1):
        activeEdge = nextEdge
        nextEdge = set()
        for pt in activeEdge:
            for n in getNeighbors(pt[0], pt[1]):
                if n not in workingSet and width > n[0] > 0 and height > n[1] > 0:  # TODO: check in bounds
                    workingSet.add(n)
                    nextEdge.add(n)
        layers.append(nextEdge)

    return layers


class AnnotatedRegion(GFF.Annotation):
    def __init__(self, GFF_annotation, renderer, xy_seq_start):
        assert isinstance(GFF_annotation, GFF.Annotation), "This isn't a proper GFF object"
        g = GFF_annotation  # short name
        super(AnnotatedRegion, self).__init__(g.chromosome, g.ID, g.source, g.feature,
                                              g.start, g.end, g.score, g.strand, g.frame,
                                              g.attributes, g.line)
        self.points = annotation_points(GFF_annotation, renderer, xy_seq_start)
        self.xy_seq_start = xy_seq_start
        self.protein_spans = []

    def exon_region_points(self):
        exon_indices = set()
        for exon in self.protein_spans:
            exon_indices.update(exon.set_of_points())
        length = len(self.points)
        exon_points = []
        for i in exon_indices:
            adjusted = i - self.start
            if 0 <= adjusted < length:
                exon_points.append(self.points[adjusted])
        return exon_points

    def add_cds_region(self, annotation_entry):
        """ :type annotation_entry: GFF.Annotation """
        self.protein_spans.append(Span(annotation_entry.start, annotation_entry.end))