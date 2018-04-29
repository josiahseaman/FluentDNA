import os
from collections import defaultdict, namedtuple
from math import log, ceil

import numpy
from DNASkittleUtils.Contigs import read_contigs, write_contigs_to_file, Contig

from DDV.DDVUtils import viridis_palette
from DDV.TileLayout import TileLayout

OligSites = namedtuple("OligSites", ['seq', 'indices'])
MIN_COPIES = 40

def find_all(a_str, sub):
    """Written by Karl Knechtel  https://stackoverflow.com/a/4665027/3067894"""
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

    #list(find_all('spam spam spam spam', 'spam')) # [0, 5, 10, 15]


def collect_divergence_snips(uniq_seq, hit, index_digits):
    snips = []
    padding = '0' + str(index_digits)
    length = len(hit.seq)
    matches = hit.indices  # tuple(find_all(ref_seq, olig))
    snips.append((hit.seq + ('-' * index_digits)) * 10)  # header with ten copies of the target olig
    for pos in matches:
        snips.append(uniq_seq[pos:pos + length])
        snips.append(format(pos, padding))
    return snips


def contig_for_each_oligsite(query_seq, oligsites, index_digits):
    """Create FASTA based on clipping out differences"""
    assert all([len(oligsites[0].seq) == len(olig.seq) for olig in oligsites]),  "All need to be the same length"

    divergence_contigs = []
    for hit in oligsites:
        divergence_snips = collect_divergence_snips(query_seq, hit, index_digits)
        if divergence_snips is not None:  # actual examples found
            divergence_contigs.append(Contig(hit.seq, ''.join(divergence_snips)))
    return divergence_contigs


def collect_likely_repeat_oligs(ref_fasta, base_width):
    ref_contigs = read_contigs(ref_fasta)
    seq = ref_contigs[0].seq
    oligs = defaultdict(lambda: 0)
    for i in numpy.random.randint(0, len(seq), len(seq) // 20):
        x = seq[ i : i + base_width]
        if '-' not in x:
            oligs[x] += 1  # increase count
    print("Found Oligs: ", len(oligs))
    threshold = 4
    candidates = [x for x in oligs if oligs[x] >= threshold]
    print(len(candidates), "Qualify for Threshold of", threshold, )
    return candidates


def collect_repeat_oligsites(seq, base_width):
    oligs = defaultdict(lambda: list())
    for i in range(len(seq) - base_width + 1):
        x = seq[i: i + base_width]
        if '-' not in x and 'N' not in x:
            oligs[x].append(i)
    passing = [OligSites(kmer, indices) for kmer, indices in oligs.items() if len(indices) >= MIN_COPIES]
    print(len(passing), "out of", len(seq), "passed")
    return passing


def visual_digit_colors(base=10):
    """Grab 0-9 digits from viridis and add them to the palette for visual indices.
    This allows you to add base 10 numbers into your sequences and they'll make some amount of sense
    even at single pixel resolution where text would be too hard to make out.
        """
    spectrum = viridis_palette()
    palette = {}
    for numeral, value in [('0', 0), ('1', 1), ('2', 2), ('3', 3), ('4', 4), ('5', 5), ('6', 6),
                           ('7', 7), ('8', 8), ('9', 9)]:
        palette[numeral] = spectrum[value * 28]  # step size of 28 gets to 9 * 28 = 252
    return palette


class DivergencePlot(TileLayout):
    """Example command:
    fluentdna.py --layout=divergence "--fasta=data\Divergence example\chrI_ce6_gapped.fa"
    --extrafastas "data\Divergence example\C_to_ce6_chrI_unique.fa" "--outname=Divergence Plot"
    --natural_colors --base_width=10"""
    def __init__(self, use_fat_headers=False, use_titles=False, sort_contigs=False, low_contrast=False,
                 base_width=10):
        self.target_oligs = []
        super().__init__(use_fat_headers, use_titles, sort_contigs, low_contrast,
                         base_width=base_width)
        self.palette.update(visual_digit_colors())


    def show_divergence(self, output_dir, output_name, ref_fasta, query_unique_fasta):
        # self.target_oligs = collect_likely_repeat_oligs(ref_fasta, self.base_width)
        ref_seq = read_contigs(ref_fasta)[0].seq
        query_seq = read_contigs(query_unique_fasta)[0].seq
        self.target_oligs = collect_repeat_oligsites(ref_seq, self.base_width)
        index_digits = ceil(log(len(query_seq)) / log(10))  # doesn't have to be base 10
        self.base_width += index_digits  # actual raster width includes the index string
        self.calculate_layout()
        divergence_contigs = contig_for_each_oligsite(query_seq, self.target_oligs, index_digits)
        divergence_file_path = os.path.join(output_dir, 'divergence.fa')
        write_contigs_to_file(divergence_file_path, divergence_contigs)
        # Display generated fasta normally
        super().process_file(divergence_file_path, output_dir, output_name)
