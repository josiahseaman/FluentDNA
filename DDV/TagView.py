from  os.path import join, basename

import math
from DNASkittleUtils.Contigs import read_contigs
from DNASkittleUtils.DDVUtils import rev_comp

from DDV.Annotations import GFF
from DDV.AnnotatedGenome import AnnotatedGenomeLayout

class TagView(AnnotatedGenomeLayout):
    def __init__(self, fasta_file, gff_file, *args, **kwargs):
        super(TagView, self).__init__(fasta_file, gff_file, *args, **kwargs)
        self.scan_width = 100
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

        # progress = 0
        # for y, line in enumerate(oligomer_profiles):
        #     if y + self.scan_width < len(oligomer_profiles):
        #         for x in range(1, self.scan_width + 1):
        #             score = pythonCorrelate(line, oligomer_profiles[y + x])
        #             match_bytes[progress] = max(0, int(score * 255))  # don't care about negative values
        #             progress += 1
        return match_bytes


    def output_tag_byte_sequence(self, output_folder, output_file_name, match_bytes):
        bytes_file = join(output_folder, output_file_name + '__100.bytes')

        with open(bytes_file, 'wb') as out:
            out.write(match_bytes)
        print("Done with line correlations")
        return bytes_file


    def render_genome(self, output_folder, output_file_name):
        # empty annotation

        # sequence
        self.contigs = read_contigs(self.fasta_file)

        # tags is viridis
        match_bytes = self.scan_reverse_complements()
        bytes_file = self.output_tag_byte_sequence(output_folder, output_file_name, match_bytes)

        ### AnnotatedGenome temporary code  ###
        # annotation_fasta = join(output_folder, basename(self.gff_filename) + '.fa')
        # chromosomes = [x.name.split()[0] for x in self.contigs]
        # lengths = [len(x.seq) for x in self.contigs]
        # create_fasta_from_annotation(self.annotation, chromosomes,
        #                              scaffold_lengths=lengths,
        #                              output_path=annotation_fasta)
        super(TagView, self).process_file(output_folder,
                          output_file_name=output_file_name,
                          fasta_files=[#annotation_fasta,
                                       self.fasta_file,
                                       bytes_file])

        pass



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







#
# def average(values, start=0, length=-1):
#     if length == -1:
#         length = len(values)
#         assert length > 0
#         return float(sum(values)) / length
#     assert length > 0
#     totalSum = 0
#     for index in range(start, start + length):
#         totalSum += values[index]
#     return float(totalSum) / length
#
#
# def pythonCorrelate(x, y):
#     n = len(x)
#     avg_x = average(x)
#     avg_y = average(y)
#     diffprod = 0.0
#     xdiff2 = 0.0
#     ydiff2 = 0.0
#     for idx in range(n):
#         xdiff = float(x[idx]) - avg_x
#         ydiff = float(y[idx]) - avg_y
#         diffprod += xdiff * ydiff
#         xdiff2 += xdiff * xdiff
#         ydiff2 += ydiff * ydiff
#
#     if xdiff2 * ydiff2 == 0.0:
#         return 0.0
#
#     return diffprod / math.sqrt(xdiff2 * ydiff2)
#


# def pearsonCorrelation(x, y):
#     """Pearson correlation coefficient between signals x and y."""
#     assert len(x) == len(y), (len(x), " vs. ", len(y))
#     n = len(x)
#     assert n > 0, "Array is empty"
#     assert isinstance(x[0], Number), x[0]
#
#     if usingCcode:
#         arrX = (ctypes.c_double * len(x))(*x)
#         arrY = (ctypes.c_double * len(y))(*y)
#
#         skittleUtils.Correlate.restype = ctypes.c_double
#         temp = skittleUtils.Correlate(arrX, arrY, n)
#         return temp
#     else:
#         return pythonCorrelate(x, y)
#
#
#
#
# BASE_DIR = r'D:\josiah\Projects\DDV\DDV'
# try:
#     if sys.platform == 'win32':
#         skittleUtils = ctypes.CDLL(os.path.join(BASE_DIR, 'libSkittleGraphUtils.dll'))
#         usingCcode = True
#         print("Optimized Windows C code for correlations found!")
#     elif 'linux' in sys.platform:
#         skittleUtils = ctypes.CDLL(os.path.join(BASE_DIR,'libSkittleGraphUtils.so.1.0.0'))
#         usingCcode = True
#         print("Optimized Linux C code for correlations found!")
#     else:
#         usingCcode = False
# except:
#     usingCcode = False
#     print("Could not find Optimized Windows C code for correlations!")
