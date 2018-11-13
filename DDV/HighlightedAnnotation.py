from itertools import chain
import traceback

import sys
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

def annotation_points(entry, renderer, start_offset):
    annotation_points = []  # TODO use unsigned shorts (max 65535) for memory
    for i in range(entry.start, entry.end):
        # important to include title and reset padding in coordinate frame
        position = renderer.relative_position(i + start_offset)
        annotation_points.append((position[0] + renderer.border_width,
                                  position[1] + renderer.border_width))
    return annotation_points


class HighlightedAnnotation(TileLayout):
    def __init__(self, gff_file, query=None, repeat_annotation=None, **kwargs):
        super(HighlightedAnnotation, self).__init__(border_width=12, **kwargs)
        self.annotation = GFF(gff_file).annotations if gff_file is not None else None
        self.query_annotation = GFF(query).annotations if query is not None else None
        self.repeat_annotation = GFF(repeat_annotation).annotations if repeat_annotation is not None else None
        self.pil_mode = 'RGBA'  # Alpha channel necessary for outline blending
        self.font_name = "ariblk.ttf"  # TODO: compatibility testing with Mac

    def process_file(self, input_file_path, output_folder, output_file_name,
                     no_webpage=False, extract_contigs=None):
        if self.annotation is not None:
            with open(input_file_path, 'r') as fasta:
                assert fasta.readline().startswith('>'), "Fasta file must start with a header '>name'"
        super(HighlightedAnnotation, self).process_file(input_file_path, output_folder, output_file_name,
                                                        no_webpage, extract_contigs)
        # nothing extra

    def draw_extras(self):
        """Drawing Annotations labels and shadow outlines"""

        positions = self.contig_struct()
        for sc_index, coordinate_frame in enumerate(positions):  # Exact match required (case sensitive)
            scaff_name = coordinate_frame["name"].split()
            if not scaff_name:
                raise ValueError('Annotation cannot proceed without a contig name in the FASTA file.  \n'
                                 'Please add a line at the beginning of your fasta file with a name '
                                 'that exactly matches the first column in your annotation. For example: '
                                 '>chrMt')
            scaff_name = scaff_name[0]
            # for one chromosome
            try:
                self.draw_extras_for_chromosome(scaff_name, coordinate_frame)
            except BaseException as e:
                print("Encountered error while rendering", scaff_name)
                traceback.print_exc()
                print(e)
                print("Continuing to next scaffold.")


    def draw_extras_for_chromosome(self, scaff_name, coordinate_frame):
        genic_color = (255, 255, 255, 46)  # faint highlighter for genic regions
        if self.repeat_annotation is not None:
            self.draw_annotation_layer(self.repeat_annotation, scaff_name, coordinate_frame, (0, 0, 0, 55),
                                       (255, 255, 255, 0), simple_entry=True)
        if self.query_annotation is not None:
            self.draw_annotation_layer(self.query_annotation, scaff_name, coordinate_frame, genic_color,
                                       (50, 50, 50, 255), shadows=True)
        if self.annotation is not None:
            self.draw_annotation_layer(self.annotation, scaff_name, coordinate_frame, genic_color,
                                       (50, 50, 50, 255))


    def draw_annotation_layer(self, annotations, scaff_name, coordinate_frame, color, label_color,
                              simple_entry=False, shadows=False):
        regions = self.find_annotated_regions(annotations, scaff_name,
                                              coordinate_frame["title_padding"], no_structure=simple_entry)
        if not len(regions):
            return  # no work to do for this scaffold

        try:
            upper_left = self.position_on_screen(coordinate_frame["xy_title_start"])
        except IndexError:
            print("Warning: Annotation layer positioning error at:", coordinate_frame["xy_title_start"],
                  file=sys.stderr)
            upper_left = [self.border_width, self.border_width]
        # relative coordinates
        width = max(set(p[0] for region in regions for p in region.points)) + 2 * self.border_width
        height = max(set(p[1] for region in regions for p in region.points)) + 2* self.border_width
        markup_image = Image.new('RGBA', (width + 1, height + 1), (0, 0, 0, 0))

        self.draw_annotation_outlines(regions, markup_image, color,
                                      simple_entry=simple_entry, shadows=shadows)

        if self.use_titles and label_color[3]:  # if the text color is transparent, don't bother
            universal_prefix = find_universal_prefix(regions)
            print("Removing Universal Prefix from annotations: '%s'" % universal_prefix)
            self.draw_annotation_labels(markup_image, regions, label_color, universal_prefix,
                                        use_suppression=simple_entry)  # labels on top
        self.image.paste(markup_image, (upper_left[0] - self.border_width,
                                        upper_left[1] - self.border_width), markup_image)


    def draw_annotation_outlines(self, regions, markup_image, color, simple_entry, shadows):
        markup_canvas = markup_image.load()
        self.draw_exons(markup_canvas, regions, color, highlight_whole_entry=True)
        if not simple_entry:
            exon_color = (255, 255, 255, 67)  # white highlighter.  This is less disruptive overall
            self.draw_exons(markup_canvas, regions, exon_color)  # double down on alpha
        if shadows:
            try:
                annotation_point_union = self.draw_big_shadow_outline(markup_image, regions, (65, 42, 80))
                self.draw_overlap_shadows(annotation_point_union, markup_image, regions, (65, 42, 80))
            except MemoryError as e:  # the global union takes a lot of memory
                print("Ran out of Memory rendering annotation shadows.  Continuing...")
                print(e)

    def draw_big_shadow_outline(self, markup_image, regions, shadow):
        print("Drawing annotation outlines")
        annotation_point_union = set()
        for region in regions:
            annotation_point_union.update(region.points)
        # desaturated purple drop shadow, decreasing opacity
        opacities = linspace(197, 10, self.border_width)
        outline_colors = [(shadow[0], shadow[1], shadow[2], int(opacity)) for opacity in opacities]
        big_shadow = outlines(annotation_point_union,
                              self.border_width, markup_image.width, markup_image.height)
        self.draw_shadow(big_shadow, markup_image.load(), outline_colors)
        return annotation_point_union

    def draw_exons(self, markup_canvas, regions, color, highlight_whole_entry=False):
        print("Drawing exons" if not highlight_whole_entry else "Drawing genic regions")

        for region in regions:
            if highlight_whole_entry:
                for point in region.points:  # highlight exons
                    blend_pixel(markup_canvas, point, color)
            else:
                for point in region.exon_region_points():  # highlight exons
                    blend_pixel(markup_canvas, point, color)

    def draw_overlap_shadows(self, annotation_point_union, markup_image, regions, shadow):
        """Find subset of genes who are completely overshadowed"""
        print("Drawing secondary shadows")
        opacities = linspace(170, 40, self.border_width // 4)
        outline_colors = [(shadow[0], shadow[1], shadow[2], int(opacity)) for opacity in opacities]
        for region in regions:
            # if annotation_point_union.issuperset(region.outline_points):
            region.outline_points = outlines(region.points,  # small outline
                                             self.border_width // 4, markup_image.width, markup_image.height)
            # remove this line if you prefer shadow intersection
            self.draw_shadow(region.outline_points, markup_image.load(), outline_colors)
            # lost_edge = region.outline_points[0].intersection(annotation_point_union)
            # if lost_edge:
            #     lost_shadows = [layer.intersection(annotation_point_union) for layer in region.outline_points]
            #     self.draw_shadow(lost_shadows, markup_canvas, outline_colors)


    def draw_shadow(self, shadow, markup_canvas, outline_colors, flat_color=False):
        for radius, layer in enumerate(shadow):
            darkness = radius
            # self.border_width - len(region.outline_points) + radius  # softer line for small features
            c = outline_colors[darkness]
            for pt in layer:
                blend_pixel(markup_canvas, pt, c)

    def find_annotated_regions(self, annotations, scaff_name, start_offset, no_structure=False):
        """:param start_offset:
           :param scaffold_name:
           :type annotations: dict(GFF.Annotation)
        """
        print("Collecting points in annotated regions of", scaff_name)
        regions = []
        genes_seen = set()
        if annotations is None:
            return regions
        if scaff_name in annotations.keys():
            for entry in annotations[scaff_name]:
                # redundancy checks for file with both mRNA and gene
                try:
                    if no_structure:  # life is simple
                        regions.append(AnnotatedRegion(entry, self, start_offset))
                    else:  # find gene/mRNA/exon/CDS hierarchy
                        if entry.feature == 'gene':
                            regions.append(AnnotatedRegion(entry, self, start_offset))
                            genes_seen.add(extract_gene_name(entry).replace('g','t'))
                        if entry.feature == 'mRNA' and extract_gene_name(entry) not in genes_seen:
                            regions.append(AnnotatedRegion(entry, self, start_offset))
                        if entry.feature == 'CDS' or entry.feature == 'exon':
                            # hopefully mRNA comes first in the file
                            # if extract_gene_name(regions[-1]) == entry.attributes['Parent']:
                            # the if was commented out because CDS [Parent] to mRNA, not gene names
                                regions[-1].add_cds_region(entry)
                except (IndexError, KeyError, TypeError, ValueError) as e:
                    print(e)
        return regions


    def draw_annotation_labels(self, markup_image, annotated_regions, label_color,
                               universal_prefix='', use_suppression=False):
        """        :type use_suppression: bool supress lines of text that would overlap and crowd
 :type annotated_regions: list(AnnotatedRegion)
        """
        print("Drawing annotation labels")
        self.fonts = {9: ImageFont.load_default()}  # clear font cache, this may be a different font
        last_unsuppressed_progress = 0
        suppression_size = 900 if use_suppression else 0
        for region in annotated_regions:
            if last_unsuppressed_progress \
                    and abs(region.start - last_unsuppressed_progress) < suppression_size:
                continue #skip
            else:
                last_unsuppressed_progress = region.start
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
                                 upper_left, vertical_label, region.strand, markup_image,
                                 label_color=label_color)
            except BaseException as e:
                print('Error while drawing label %s' % extract_gene_name(region), e)

    def handle_multi_column_annotations(self, region, left, right, top, bottom):
        multi_column = abs(right - left) > self.base_width
        if multi_column:  # pick the biggest column to contain the label, ignore others
            median_point = len(region.points) // 2 + min(region.start, region.end)
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
                    strand, canvas, horizontal_centering=False, center_vertical=False, chop_text=True,
                    label_color=(50, 50, 50, 255)):
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
        if font_size >= 14:
            alpha = 235 / 255
            if font_size > 30:
                alpha = 200 / 255
            label_color = (label_color[0], label_color[1], label_color[2], int(label_color[3] * alpha))
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
        del txt


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
                if n not in workingSet \
                        and width > n[0] > 0 and height > n[1] > 0:  # TODO: check in bounds
                    workingSet.add(n)
                    nextEdge.add(n)
        layers.append(nextEdge)

    return layers


class AnnotatedRegion(GFF.Annotation):
    def __init__(self, GFF_annotation, renderer, start_offset):
        assert isinstance(GFF_annotation, GFF.Annotation), "This isn't a proper GFF object"
        g = GFF_annotation  # short name
        super(AnnotatedRegion, self).__init__(g.chromosome, g.ID, g.source, g.feature,
                                              g.start, g.end, g.score, g.strand, g.frame,
                                              g.attributes, g.line)
        self.points = annotation_points(GFF_annotation, renderer, start_offset)
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