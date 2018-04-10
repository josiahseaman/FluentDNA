import os

import numpy
from DNASkittleUtils.Contigs import read_contigs, write_contigs_to_file, Contig

from DDV.TileLayout import TileLayout


def find_all(a_str, sub):
    """Written by Karl Knechtel  https://stackoverflow.com/a/4665027/3067894"""
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

    #list(find_all('spam spam spam spam', 'spam')) # [0, 5, 10, 15]


def collect_divergence_snips(ref_seq, uniq_seq, olig):
    snips = []
    length = len(olig)
    matches = tuple(find_all(ref_seq, olig))
    if len(matches) > 40:  # actual examples found
        snips.append(olig*10)  # header with ten copies of the target olig
        for pos in matches:
            snips.append(uniq_seq[pos:pos + length])
        return snips
    else:
        return None


def contig_for_each_olig(ref_fasta, query_unique_fasta, oligs):
    """Create FASTA based on clipping out differences"""
    assert all([len(oligs[0]) == len(olig) for olig in oligs]),  "All need to be the same length"
    ref_contigs = read_contigs(ref_fasta)
    uniq_contigs = read_contigs(query_unique_fasta)
    divergence_contigs = []
    for olig in oligs:
        divergence_snips = collect_divergence_snips(ref_contigs[0].seq, uniq_contigs[0].seq, olig)
        if divergence_snips is not None:  # actual examples found
            divergence_contigs.append(Contig(olig, ''.join(divergence_snips)))
    return divergence_contigs




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

    def show_divergence(self, output_dir, output_name, ref_fasta, query_unique_fasta):
        self.collect_1000_oligs(ref_fasta)
        divergence_contigs = contig_for_each_olig(ref_fasta, query_unique_fasta, self.target_oligs)
        divergence_file_path = os.path.join(output_dir, 'divergence.fa')
        write_contigs_to_file(divergence_file_path, divergence_contigs)
        # Display generated fasta normally
        super().process_file(divergence_file_path, output_dir, output_name)

    def collect_1000_oligs(self, ref_fasta):
        ref_contigs = read_contigs(ref_fasta)
        seq = ref_contigs[0].seq
        oligs = set()
        for i in numpy.random.randint(0, len(seq), 100000):
            x = seq[ i : i + self.base_width]
            if '-' not in x:
                oligs.add(x)
        print("Found Oligs: ", len(oligs))
        self.target_oligs = list(oligs)