import traceback
from datetime import datetime
from  os.path import join, basename

import math
from DNASkittleUtils.Contigs import read_contigs
from PIL import Image

from DDV.Annotations import GFF
from DDV.TileLayout import TileLayout, hex_to_rgb
from collections import namedtuple
Point = namedtuple('Point', ['x', 'y'])

class OutlinedAnnotation(TileLayout):
    def __init__(self, fasta_file, gff_file, **kwargs):
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

        self.draw_annotation_outlines(markup_canvas)
        self.draw_annotation_labels(markup_canvas)
        self.image = Image.alpha_composite(self.image, markup_image)

    def draw_annotation_outlines(self, markup_canvas):
        regions = self.find_annotated_regions()
        print("Drawing annotation outlines")
        outline_colors = [(58, 20, 84, 197),  # desaturated purple drop shadow, decreasing opacity
                          (58, 20, 84, 162),
                          (58, 20, 84, 128),
                          (58, 20, 84, 84),
                          (58, 20, 84, 39),
                          (58, 20, 84, 15)]
        for region in regions:
            for radius, layer in enumerate(region.outline_points):
                c = outline_colors[radius]
                for pt in layer:
                    if markup_canvas[pt[0], pt[1]][3] == 0:  # nothing drawn
                        markup_canvas[pt[0], pt[1]] = c
                    else:
                        remaining_light = 1.0 - (markup_canvas[pt[0], pt[1]][3] / 256)
                        combined_alpha = 256 - int(remaining_light * (256 - c[3]) )
                        markup_canvas[pt[0], pt[1]] = (c[0], c[1], c[2], combined_alpha)

    def find_annotated_regions(self):
        print("Collecting points in annotated regions")
        positions = self.contig_struct()
        regions = []
        for sc_index, coordinate_frame in enumerate(positions):  # Exact match required (case sensitive)
            scaff_name = coordinate_frame["name"].split()[0]
            if scaff_name in self.annotation.annotations.keys():
                for entry in self.annotation.annotations[scaff_name]:
                    if entry.feature == 'gene':
                        annotation_points = []  # this became too complex for a list comprehension
                        for i in range(entry.start, entry.end):
                            # important to include title and reset padding in coordinate frame
                            progress = i + coordinate_frame["xy_seq_start"]
                            annotation_points.append(self.position_on_screen(progress))
                        regions.append(AnnotatedRegion(entry, annotation_points))
        return regions


    def draw_annotation_labels(self, markup_canvas):
        print("Drawing annotation labels")


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
        self.points = set(annotation_points)
        self.outline_points = outlines(annotation_points, 6)
