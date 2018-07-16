import traceback
from datetime import datetime
from  os.path import join, basename

import math
from DNASkittleUtils.Contigs import read_contigs
from PIL import Image

from DDV.Annotations import GFF
from DDV.Span import Span
from DDV.TileLayout import TileLayout, hex_to_rgb
from collections import namedtuple
Point = namedtuple('Point', ['x', 'y'])


def blend_pixel(markup_canvas, pt, c):
    if markup_canvas[pt[0], pt[1]][3] == 0:  # nothing drawn
        markup_canvas[pt[0], pt[1]] = c
    else:
        remaining_light = 1.0 - (markup_canvas[pt[0], pt[1]][3] / 256)
        combined_alpha = 256 - int(remaining_light * (256 - c[3]) )
        markup_canvas[pt[0], pt[1]] = (c[0], c[1], c[2], combined_alpha)


class OutlinedAnnotation(TileLayout):
    def __init__(self, fasta_file, gff_file, **kwargs):
        kwargs['font_name'] ="ariblk.ttf"
        super(OutlinedAnnotation, self).__init__(**kwargs)
        self.fasta_file = fasta_file
        self.gff_filename = gff_file
        self.annotation = GFF(self.gff_filename)
        self.pil_mode = 'RGBA'  # Alpha channel necessary for outline blending


    def process_file(self, input_file_path, output_folder, output_file_name):
        super(OutlinedAnnotation, self).process_file(input_file_path, output_folder, output_file_name)
        # nothing extra

    def draw_titles(self):
        super(OutlinedAnnotation, self).draw_titles()
        markup_image = Image.new('RGBA', (self.image.width, self.image.height), (0,0,0,0))
        markup_canvas = markup_image.load()

        annotated_regions = self.draw_annotation_outlines(markup_canvas)
        self.draw_annotation_labels(markup_image, annotated_regions)
        self.image = Image.alpha_composite(self.image, markup_image)

    def draw_annotation_outlines(self, markup_canvas):
        regions = self.find_annotated_regions()
        print("Drawing annotation outlines")
        outline_colors = [(58, 20, 84, 197),  # desaturated purple drop shadow, decreasing opacity
                          (58, 20, 84, 162),
                          (58, 20, 84, 128),
                          (58, 20, 84, 84),
                          (58, 20, 84, 49),
                          (58, 20, 84, 15)]
        exon_color = (255,255,255,135)  # white highlighter.  This is less disruptive overall
        for region in regions:
            for radius, layer in enumerate(region.outline_points):
                darkness = 6 - len(region.outline_points) + radius  # softer line for small features
                c = outline_colors[darkness]
                for pt in layer:
                    blend_pixel(markup_canvas, pt, c)

            for point in region.dark_region_points():
                blend_pixel(markup_canvas, point, exon_color)
        return regions

    def find_annotated_regions(self):
        print("Collecting points in annotated regions")
        positions = self.contig_struct()
        regions = []
        for sc_index, coordinate_frame in enumerate(positions):  # Exact match required (case sensitive)
            scaff_name = coordinate_frame["name"].split()[0]
            if scaff_name in self.annotation.annotations.keys():
                for entry in self.annotation.annotations[scaff_name]:
                    if entry.feature == 'mRNA':
                        annotation_points = []  # this became too complex for a list comprehension
                        for i in range(entry.start, entry.end):
                            # important to include title and reset padding in coordinate frame
                            progress = i + coordinate_frame["xy_seq_start"]
                            annotation_points.append(self.position_on_screen(progress))
                        regions.append(AnnotatedRegion(entry, annotation_points))
                    if entry.feature == 'CDS':
                        # hopefully mRNA comes first in the file
                        if regions[-1].attributes['Name'] == entry.attributes['Parent']:
                            regions[-1].add_cds_region(entry)

        return regions


    def draw_annotation_labels(self, markup_image, annotated_regions):
        """ :type annotated_regions: list(AnnotatedRegion) """
        print("Drawing annotation labels")
        for region in annotated_regions:
            xs = [pt[0] for pt in region.points]
            left_most, right_most = min(xs), max(xs)
            ys = [pt[1] for pt in region.points]
            top, bottom = min(ys), max(ys)
            multi_column = abs(right_most - left_most) > self.base_width
            if multi_column:
                median_point = region.points[len(region.points) // 2]
            width = self.base_width
            height = len(region.points) // self.base_width
            vertical_label = height > width
            rotation_direction = -90 if region.strand == '-' else 90

            upper_left = (left_most, top)
            bottom_right = (right_most, bottom)
            # width, height = bottom_right[0] - upper_left[0], bottom_right[1] - upper_left[1]

            title_width = 18

            # Title orientation and size
            if vertical_label:
                width, height = height, width  # swap
            font_size = max(9, int((width * .03222) + 6))  # found eq with two reference points

            self.write_title(region.attributes["Name"], width, height, font_size, 1, title_width,
                             upper_left, vertical_label, markup_image)


def getNeighbors(pt):
    return {(pt[0] + 1, pt[1]), (pt[0] - 1, pt[1]), (pt[0], pt[1] + 1), (pt[0], pt[1] - 1)}

def allNeighbors(pt):
    return getNeighbors(pt).union({(pt[0] + 1, pt[1] + 1), (pt[0] - 1, pt[1] - 1),
                                   (pt[0] - 1, pt[1] + 1), (pt[0] + 1, pt[1] - 1)})

def outlines(annotation_points, radius, square_corners=False):
    workingSet = set()
    nextEdge = set()
    workingSet.update(annotation_points)
    nextEdge.update(annotation_points)
    layers = []
    for iterationStep in range(radius, 0,  -1):
        activeEdge = nextEdge
        nextEdge = set()

        for block in activeEdge:
            neighbors = allNeighbors(block) if square_corners else getNeighbors(block)
            for n in neighbors:
                if n not in workingSet and n[0] > 0 and n[1] > 0:  # TODO: check in bounds
                    workingSet.add(n)
                    nextEdge.add(n)
        layers.append(nextEdge)
    return layers


class AnnotatedRegion(GFF.Annotation):
    def __init__(self, GFF_annotation, annotation_points):
        assert isinstance(GFF_annotation, GFF.Annotation), "This isn't a proper GFF object"
        g = GFF_annotation  # short name
        super(AnnotatedRegion, self).__init__(g.chromosome, g.ID, g.source, g.feature,
                                              g.start, g.end, g.score, g.strand, g.frame,
                                              g.attributes, g.line)
        self.points = list(annotation_points)
        radius = 6 if self.feature == 'mRNA' else 3
        self.outline_points = outlines(annotation_points, radius)
        self.protein_spans = []

    def dark_region_points(self):
        exon_indices = []
        for i in range(self.start, self.end):
            if any([i in exon for exon in self.protein_spans]):
                exon_indices.append(i)
        exon_points = [self.points[i - self.start] for i in exon_indices]
        # introns = [i for span in self.non_protein_spans for i in range(span.begin, span.end)]
        return exon_points

    def add_cds_region(self, annotation_entry):
        """ :type annotation_entry: GFF.Annotation """
        self.protein_spans.append(Span(annotation_entry.start, annotation_entry.end))