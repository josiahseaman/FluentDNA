from  os.path import join, basename, exists

import math
import numpy as np
from DNASkittleUtils.Contigs import read_contigs
from DNASkittleUtils.DDVUtils import rev_comp

from DDV.Annotations import GFF, create_fasta_from_annotation
from DDV.AnnotatedGenome import AnnotatedGenomeLayout

class TagView(AnnotatedGenomeLayout):
    def __init__(self, fasta_file, ref_annotation, **kwargs):
        super(AnnotatedGenomeLayout, self).__init__(n_genomes=3, **kwargs)  # skipping parent
        self.fasta_file = fasta_file
        self.gff_filename = ref_annotation
        # Important: attribute_sep changed from AnnotatedGenomeLayout
        self.annotation = GFF(self.gff_filename, attribute_sep=' ')
        self.scan_width = self.base_width
        self.oligomer_size = 9

    def scan_reverse_complements(self):
        assert len(self.contigs) == 1, "Read one sequence first"
        seq = self.contigs[0].seq
        n_lines = (len(seq) // self.scan_width)
        match_bytes = bytearray(n_lines * self.scan_width)
        # calculate oligomer profiles
        print("Calculating line correlations")

        observedMax = max(self.scan_width / 6.0, 5.0)  # defined by observation of histograms
        observedOligsPerLine = setOfObservedOligs(seq, self.scan_width, self.oligomer_size)
        observedRevCompOligs = reverseComplementSet(observedOligsPerLine)
        for y in range(min(len(observedOligsPerLine), n_lines)):
            for x in range(1, min(len(observedOligsPerLine) - y - 1, self.scan_width)):
                if x + y < len(observedOligsPerLine):  # account for second to last chunk
                    matches = observedOligsPerLine[y].intersection(observedRevCompOligs[y + x])
                    if matches:
                        match_bytes[y * self.scan_width + x] = int(min(len(matches) / observedMax * 255, 255))

        return match_bytes


    def output_tag_byte_sequence(self, bytes_file, match_bytes):

        with open(bytes_file, 'wb') as out:
            out.write(match_bytes)
        print("Done with line correlations")
        return bytes_file


    def render_genome(self, output_folder, output_file_name):
        # empty annotation

        # sequence
        self.contigs = read_contigs(self.fasta_file)

        # tags is viridis
        bytes_file = join(output_folder, output_file_name + '__%i.bytes' % self.scan_width)
        if not exists(bytes_file):
            match_bytes = self.scan_reverse_complements()
            bytes_file = self.output_tag_byte_sequence(bytes_file, match_bytes)

        ### AnnotatedGenome temporary code  ###
        annotation_fasta = join(output_folder, basename(self.gff_filename) + '.fa')
        chromosomes = [x.name.split()[0] for x in self.contigs]
        lengths = [len(x.seq) for x in self.contigs]
        create_fasta_from_annotation(self.annotation, chromosomes,
                                     scaffold_lengths=lengths,
                                     output_path=annotation_fasta,
                     # TODO: currently default because all entries are "exon" in sample file
                                     features=None)
        super(TagView, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[annotation_fasta,
                                       self.fasta_file,
                                       bytes_file])

    def name_for_annotation_entry(self, entry):
        name = entry.attributes['gene_id']  # based on the RepeatMasker GTF that I have, ' ' separator
        return name

def hasDepth(listLike):
    try:
        if len(listLike) > 0 and not isinstance(listLike, (str, dict, tuple, type(u"unicode string"))) and hasattr(
                listLike[0], "__getitem__"):
            return True
        else:
            return False
    except:
        return False

def chunkUpList(seq, chunkSize, overlap=0):
    if hasDepth(seq):
        return map(lambda x: chunkUpList(x, chunkSize, overlap), seq)
    if chunkSize == 0:
        return seq
    height = int(math.ceil(len(seq) / float(chunkSize)))
    #    if height == 0: return []
    resultVector = [seq[chunk * chunkSize: (chunk + 1) * chunkSize + overlap] for chunk in range(height)]
    return resultVector



def reverseComplementSet(observedOligsPerLine):
    """
    :param observedOligsPerLine:
    :rtype: list of set of strings
    :return: set of reverse complements of input set
    """
    lines = []
    for oligs in observedOligsPerLine:
        lines.append({rev_comp(word) for word in oligs})
    return lines


def observedOligs(seq, oligomerSize):
    oligs = set()
    for endIndex in range(oligomerSize, len(seq) + 1, 1):
        window = seq[endIndex - oligomerSize: endIndex]
        oligs.add(window)
    return oligs


def setOfObservedOligs(seq, lineSize, oligomerSize):
    # TODO: refactor and remove duplication from OligomerUsage.countOligomers()
    """
    :return: list of set of strings
    :rtype: list
    """
    overlap = oligomerSize - 1
    # chunk sequence by display line #we can't do this simply by line because of the overhang of oligState.oligState
    lines = chunkUpList(seq, lineSize, overlap)

    oligsByLine = []
    for line in lines:
        oligsByLine.append(observedOligs(line, oligomerSize))

    return oligsByLine





