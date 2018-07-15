import traceback
from datetime import datetime
from  os.path import join, basename

import math
from DNASkittleUtils.Contigs import read_contigs

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

    def process_file(self, input_file_path, output_folder, output_file_name):
        super(OutlinedAnnotation, self).process_file(input_file_path, output_folder, output_file_name)
        # nothing extra

    def draw_titles(self):
        super(OutlinedAnnotation, self).draw_titles()
        self.draw_annotation_outlines()
        self.draw_annotation_labels()

    def draw_annotation_outlines(self):
        regions = self.find_annotated_regions()
        print("Drawing annotation outlines")
        outline_colors = [hex_to_rgb(h) for h in ('#775385', "825d8b", 'a37eab', 'ba92b5', '#c08eb8') ]  # purple desaturated
        for region in regions:
            for radius, layer in enumerate(region.outline_points):
                c = outline_colors[radius]
                for pt in layer:
                    self.pixels[pt.x, pt.y] = c

    def find_annotated_regions(self):
        print("Collecting points in annotated regions")
        positions = self.contig_struct()
        regions = []
        for sc_index, coordinate_frame in enumerate(positions):  # Exact match required (case sensitive)
            scaff_name = coordinate_frame["name"].split()[0]
            if scaff_name in self.annotation.annotations.keys():
                for entry in self.annotation.annotations[scaff_name]:
                    annotation_points = []  # this became too complex for a list comprehension
                    for i in range(entry.start):
                        # important to include title and reset padding in coordinate frame
                        progress = i + coordinate_frame["xy_seq_start"]
                        annotation_points.append(Point(*self.position_on_screen(progress)))
                    regions.append(AnnotatedRegion(entry, annotation_points))
        return regions


    def draw_annotation_labels(self):
        print("Drawing annotation labels")


def getNeighbors(pt):
    return {Point(pt.x + 1, pt.y), Point(pt.x - 1, pt.y), Point(pt.x, pt.y + 1), Point(pt.x, pt.y - 1)}

def allNeighbors(pt):
    return getNeighbors(pt).union({Point(pt.x + 1, pt.y + 1), Point(pt.x - 1, pt.y - 1),
                                   Point(pt.x - 1, pt.y + 1), Point(pt.x + 1, pt.y - 1)})

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
                if n not in workingSet and n.x > 0 and n.y > 0:  # TODO: check in bounds
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
        self.outline_points = outlines(annotation_points, 1)
