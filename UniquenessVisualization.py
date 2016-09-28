from ChainFiles import fetch_all_chains
from ChainParser import ChainParser, pluck_contig


def remove_from_range(original, remove_this):
    """Original is a range defined by (start, end).  Remove a middle range 'remove_this'
    with a (start, end) and you end up with a pair of two ranges on either side of the removal.
    Special casing for the removal overlapping the beginning or end."""
    first = (original[0], remove_this[0])
    second = (remove_this[1], original[1])

    if remove_this[0] <= original[0] and remove_this[1] >= original[1]:  # delete the whole thing
        return None, None
    if remove_this[0] <= original[0] < remove_this[1]:  # overlaps start
        return None, second
    if remove_this[1] >= original[1] > remove_this[0]:  # overlaps ending
        return first, None

    return first, second  # happy path


class UniquenessViz(ChainParser):
    def __init__(self, *args, **kwargs):
        super(UniquenessViz, self).__init__(*args, **kwargs)
        self.uncovered_areas = []  # Absolute coordinates.  highly mutable: better as a blist


    def find_zero_coverage_areas(self):
        """Start with whole chromosome, subtract coverage from there"""
        self.uncovered_areas = [(0, len(self.ref_sequence))]  # TODO: zero indexed?
        master = fetch_all_chains('chr20', 'chr20', '+', self.chain_list)[0]
        ref_pointer = master.tStart  # where we are in the reference  Todo: or is it 1?
        scrutiny_index = 0
        for entry in master.entries:
            first, second = remove_from_range(self.uncovered_areas.pop(scrutiny_index), (ref_pointer, ref_pointer + entry.size))
            while second is None:  # eating the tail, this could be multiple ranges before the end is satisfied
                if first is not None:  # leaving behind a beginning, advance scrutiny_index by 1
                    self.uncovered_areas.insert(scrutiny_index, first)
                    scrutiny_index += 1  # insert first, but the next thing to be checked is whether it affects the next area as well
                else:
                    pass  # if both are None, then we've just deleted the uncovered entry
                # check again on the next one
                first, second = remove_from_range(self.uncovered_areas.pop(scrutiny_index), (ref_pointer, ref_pointer + entry.size))

            if first is None:
                self.uncovered_areas.insert(scrutiny_index, second)
            elif None not in [first, second]:  # neither are None
                self.uncovered_areas.insert(scrutiny_index, first)
                scrutiny_index += 1  # since chain entries don't overlap themselves, we're always mutating the next entry
                self.uncovered_areas.insert(scrutiny_index, second)

            ref_pointer += entry.size + entry.gap_query


    def write_zero_coverage_areas(self, ref_name):
        uniq_collection = []  # of strings
        ref_unique_name = ref_name[:-10] + '_unique.fa'
        self.find_zero_coverage_areas()
        for region in self.uncovered_areas:
            if region[1] - region[0] > 20:
                uniq_collection.append(self.ref_sequence[region[0]: region[1]])

        with open(ref_unique_name, 'w') as ref_filestream:
            self.write_fasta_lines(ref_filestream, ''.join(uniq_collection))
        return ref_unique_name


    def read_query_seq_to_memory(self, query_chr, query_source):
        pass  # we don't actually need the query sequence, only the chain file


    def main(self, chromosome_name):
        import DDV
        names, ignored, ref_chr = self.setup_for_reference_chromosome(chromosome_name)
        self.ref_sequence = pluck_contig(ref_chr, names['ref'])  # only need the reference chromosome read, skip the others
        # actual work
        names['ref_unique'] = self.write_zero_coverage_areas(names['ref'])

        folder_name = self.output_folder_prefix + ref_chr
        source_path = '.\\bin\\Release\\output\\dnadata\\'
        if True:  #self.trial_run:  # these files are never used in the viz
            del names['ref']
            del names['query']
        self.move_fasta_source_to_destination(names, folder_name, source_path)
        DDV.DDV_main(['DDV', names['ref_unique'], source_path, folder_name])



def do_chromosome(chr):
    parser = UniquenessViz(chain_name='hg38ToPanTro4.over.chain',
                           first_source='hg38_chr20.fa',  # 'HongKong\\hg38.fa',
                           second_source='panTro4.fa',  # 'panTro4_chr20.fa',
                           output_folder_prefix='Hg38_unique20_vs_panTro4_',
                           trial_run=False)
    parser.main(chr)


if __name__ == '__main__':
    chromosomes = [('chr20', 'chr20')]  # 'chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr22 chrX'.split()
    # TODO: handle Chr2A and Chr2B separately
    for chr in chromosomes:
        do_chromosome(chr)

    # import multiprocessing
    # workers = multiprocessing.Pool(6)  # number of simultaneous processes.  Watch your RAM usage
    # workers.map(do_chromosome, chromosomes)

