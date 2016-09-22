import os
import shutil
from array import array
from collections import defaultdict
from datetime import datetime

from ChainFiles import chain_file_to_list, fetch_all_chains


def chunks(seq, size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def pluck_contig(chromosome_name, genome_source):
    """Scan through a genome fasta file looking for a matching contig name.  When it find it, find_contig collects
    the sequence and returns it as a string with no cruft."""
    chromosome_name = '>' + chromosome_name
    contig_header = ''
    seq_collection = []
    printing = False
    with open(genome_source, 'r') as genome:
        for line in genome.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                # headers.append(line)
                if line == chromosome_name:
                    printing = True
                    print("Found", line)
                    contig_header = line
                elif printing:
                    break  # we've collected all sequence and reached the beginning of the next contig
            elif printing:  # This MUST come after the check for a '>'
                seq_collection.append(line.upper())  # always upper case so equality checks work
    assert len(seq_collection), "Contig not found." + chromosome_name  # File contained these contigs:\n" + '\n'.join(headers)
    return ''.join(seq_collection), contig_header


def first_word(string):
    import re
    if '\\' in string:
        string = string[string.rindex('\\') + 1:]
    return re.split('[\W_]+', string)[0]


def rev_comp(plus_strand):
    comp = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N', 'X': 'X'}
    return ''.join([comp[a] for a in reversed(plus_strand)])


class ChainParser:
    def __init__(self, chain_name, second_source, first_source, output_folder_prefix, trial_run=False, swap_columns=False):
        self.width_remaining = defaultdict(lambda: 70)
        self.ref_source = first_source  # example hg38ToPanTro4.chain  hg38 is the reference, PanTro4 is the query (has strand flips)
        self.query_source = second_source
        self.output_folder_prefix = output_folder_prefix
        self.query_contigs = dict()
        self.trial_run = trial_run
        self.swap_columns = swap_columns
        self.query_sequence = ''
        self.ref_sequence = ''
        self.query_seq_gapped = array('u', '')
        self.ref_seq_gapped = array('u', '')
        self.relevant_chains = []  # list of contigs that contribute to ref_chr in chain entries
        self.stored_rev_comps = {}

        self.chain_list = chain_file_to_list(chain_name)

    def read_contigs(self, input_file_path):
        print("Reading contigs... ", input_file_path)
        start_time = datetime.now()
        self.query_contigs = {}
        current_name = ""
        seq_collection = []

        # Pre-read generates an array of contigs with labels and sequences
        with open(input_file_path, 'r') as streamFASTAFile:
            for read in streamFASTAFile.read().splitlines():
                if read == "":
                    continue
                if read[0] == ">":
                    # If we have sequence gathered and we run into a second (or more) block
                    if len(seq_collection) > 0:
                        sequence = "".join(seq_collection)
                        seq_collection = []  # clear
                        self.query_contigs[current_name] = sequence
                    current_name = read[1:].strip()  # remove >
                else:
                    # collects the sequence to be stored in the contig, constant time performance don't concat strings!
                    seq_collection.append(read.upper())

        # add the last contig to the list
        sequence = "".join(seq_collection)
        self.query_contigs[current_name] = sequence
        print("Read FASTA Contigs in:", datetime.now() - start_time)


    def read_seq_to_memory(self, query_chr, ref_chr, query_source, ref_source, query_file_name, ref_file_name):
        self.read_contigs(query_source)
        self.query_sequence = self.query_contigs[query_chr]

        # if not self.trial_run:
        #     self.write_fasta_lines(query_file_name, self.query_sequence)
        #     self.write_fasta_lines(ref_file_name, self.ref_sequence)


    def write_fasta_lines(self, filestream, seq):
        if isinstance(filestream, str):  # I'm actually given a file name and have to open it myself
            file_name = filestream
            with open(file_name, 'w') as filestream:
                self.do_write(filestream, 70, seq)
                print("Wrote", file_name, len(self.ref_sequence))
        else:
            remaining = self.width_remaining[filestream]
            self.do_write(filestream, remaining, seq)

    def do_write(self, filestream, remaining, seq):
        try:
            left_over_size = 70
            remainder = seq[:remaining]
            seq = seq[remaining:]
            filestream.write(remainder + '\n')
            for line in chunks(seq, 70):
                if len(line) == 70:
                    filestream.write(line + '\n')
                else:
                    left_over_size = 70 - len(line)
                    filestream.write(line)  # no newline
            self.width_remaining[filestream] = left_over_size
        except Exception as e:
            print(e)

    def mash_fasta_and_chain_together(self, chain, is_master_alignment=False):
        ref_pointer, query_pointer = self.setup_chain_start(chain, is_master_alignment)

        if not is_master_alignment:
            self.do_translocation_housework(chain, query_pointer, ref_pointer)
        self.process_chain_body(chain, query_pointer, ref_pointer, is_master_alignment)


    def process_chain_body(self, chain, query_pointer, ref_pointer, is_master_alignment):
        for entry in chain.entries:  # ChainEntry
            size, gap_query, gap_ref = entry.size, entry.gap_query, entry.gap_ref

            # Debugging code
            if self.trial_run and len(self.ref_seq_gapped) > 1000000:  # 9500000  is_master_alignment and
                break

            query_seq_absolute = self.query_sequence[query_pointer: query_pointer + size + gap_ref]
            query_snippet = query_seq_absolute + 'X' * gap_query
            query_pointer += size + gap_ref  # alignable and unalignable block concatenated together
            self.query_seq_gapped.extend(query_snippet)

            ref_snippet = self.ref_sequence[ref_pointer: ref_pointer + size] + 'X' * gap_ref
            ref_snippet += self.ref_sequence[ref_pointer + size: ref_pointer + size + gap_query]
            ref_pointer += size + gap_query  # two blocks of sequence separated by gap
            self.ref_seq_gapped.extend(ref_snippet)

            if len(ref_snippet) != len(query_snippet):
                print(len(query_snippet), len(ref_snippet), "You should be outputting equal length strings til the end")

        # if is_master_alignment:  # last one: print out all remaining sequence
        #     # query_gapped_strand.extend(self.query_sequence[query_pointer:])
        #     # self.ref_seq_gapped.extend(self.ref_sequence[ref_pointer:])
        # else:


    def setup_chain_start(self, chain, is_master_alignment):
        # convert to int except for chr names and strands
        ref_pointer, query_pointer = chain.tStart, chain.qStart
        if chain.tEnd - chain.tStart > 100 * 1000:
            print('>>>>', chain)
        if is_master_alignment:  # include the unaligned beginning of the sequence
            longer_gap = max(query_pointer, ref_pointer)
            self.query_seq_gapped.extend(self.query_sequence[:query_pointer] + 'X' * (longer_gap - query_pointer))  # one of these gaps will be 0
            self.ref_seq_gapped.extend(self.ref_sequence[:ref_pointer] + 'X' * (longer_gap - ref_pointer))
        return ref_pointer, query_pointer


    def do_translocation_housework(self, chain, query_pointer, ref_pointer):
        minus_strand = chain.qStrand != '+'
        filler = 'G' if not minus_strand else 'C'
        self.query_seq_gapped.extend(filler * (500 + (100 - len(self.query_seq_gapped) % 100)))  # visual separators
        self.ref_seq_gapped.extend(filler * (500 + (100 - len(self.ref_seq_gapped) % 100)))
        # label, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = header.split()
        # self.query_seq_gapped.extend('_'.join(['\n>' + tName, qStrand, qName]) + '\n')  # visual separators
        # self.ref_seq_gapped.extend('_'.join(['\n>' + qName, qStrand, tName]) + '\n')

        # delete the ungapped query sequence
        # 	delete the query sequence that doesn't match to anything based on the original start, stop, size,
        # 	replace query with X's, redundant gaps will be closed later
        # 		there probably shouldn't be any gap at all between where the gap starts and where it gets filled in by the new chain
        # 	what if that overlaps to reference sequence and not just N's?

        # delete the target reference region
        # 	target reference will be filled in with gapped version of reference (might be slightly longer)
        # 	compensate start position for all previous gaps in the reference
        # 	delete parallel query region (hopefully filled with N's)

        # insert the gapped versions on both sides
        pass

    def write_gapped_fasta(self, query, reference):
        query_gap_name = os.path.splitext(query)[0] + '_gapped.fa'
        ref_gap_name = os.path.splitext(reference)[0] + '_gapped.fa'
        with open(query_gap_name, 'w') as query_file:
            self.write_fasta_lines(query_file, ''.join(self.query_seq_gapped))
        with open(ref_gap_name, 'w') as ref_file:
            self.write_fasta_lines(ref_file, ''.join(self.ref_seq_gapped))
        return query_gap_name, ref_gap_name


    def print_only_unique(self, query_gapped_name, ref_gapped_name):
        query_uniq_array = array('u', self.query_seq_gapped)
        que_uniq_array = array('u', self.ref_seq_gapped)
        print("Done allocating unique array")
        shortest_sequence = min(len(self.query_seq_gapped), len(self.ref_seq_gapped))
        # scanning_past_header = False
        for i in range(shortest_sequence):
            # if self.query_seq_gapped[i] == '\n':
            #     scanning_past_header = False  # stop scanning once you hit the terminating newline
            #     continue
            # if self.query_seq_gapped[i] in '>;':  # ; is for comments
            #     scanning_past_header = True
            # if scanning_past_header:  # query_uniq_array is already initialized to contain header characters
            #     continue
            # only overlapping section
            if self.query_seq_gapped[i] == self.ref_seq_gapped[i]:
                query_uniq_array[i] = 'X'
                que_uniq_array[i] = 'X'
            else:  # Not equal
                if self.query_seq_gapped[i] == 'N':
                    query_uniq_array[i] = 'X'
                if self.ref_seq_gapped[i] == 'N':
                    que_uniq_array[i] = 'X'

        # Just to be thorough: prints aligned section (shortest_sequence) plus any dangling end sequence
        query_unique_name = query_gapped_name[:-10] + '_unique.fa'
        ref_unique_name = ref_gapped_name[:-10] + '_unique.fa'
        with open(query_unique_name, 'w') as query_filestream:
            self.write_fasta_lines(query_filestream, ''.join(query_uniq_array))
        with open(ref_unique_name, 'w') as ref_filestream:
            self.write_fasta_lines(ref_filestream, ''.join(que_uniq_array))

        return query_unique_name, ref_unique_name


    @staticmethod
    def move_fasta_source_to_destination(fasta, folder_name, source_path):
        destination_folder = os.path.join(source_path, folder_name)
        if os.path.exists(destination_folder):
            shutil.rmtree(destination_folder, ignore_errors=True)  # Make sure we can overwrite the contents
        os.makedirs(destination_folder, exist_ok=False)
        for key in fasta:
            shutil.move(fasta[key], destination_folder)
            fasta[key] = os.path.join(destination_folder, fasta[key])


    def do_all_relevant_chains(self):
        """In panTro4ToHg38.over.chain there are ZERO chains that have a negative strand on the reference 'tStrand'.
        I think it's a rule that you always flip the query strand instead."""
        previous = None
        for chain in self.relevant_chains:
            if chain.qStrand == '-':  # need to load rev_comp
                if chain.qName not in self.stored_rev_comps:
                    if len(self.query_contigs[chain.qName]) > 1000000:
                        print("Reversing", chain.qName, len(self.query_contigs[chain.qName]) // 1000000, 'Mbp')
                    self.stored_rev_comps[chain.qName] = rev_comp(self.query_contigs[chain.qName])  # caching for performance
                    # TODO: replace rev_comp with a lazy generator that mimics having the whole string using [x] and [y:x]
                self.query_sequence = self.stored_rev_comps[chain.qName]
            else:
                self.query_sequence = self.query_contigs[chain.qName]
            self.mash_fasta_and_chain_together(chain, previous is None)
            previous = chain


    def do_all_chain_matches(self, ref_chr, query_chr, query_strand):
        """In panTro4ToHg38.over.chain there are ZERO chains that have a negative strand on the reference 'tStrand'.
        I think it's a rule that you always flip the query strand instead."""
        translocations = fetch_all_chains(ref_chr, query_chr, query_strand, self.chain_list)
        if len(translocations):
            if query_strand == '-':
                self.query_sequence = rev_comp(self.query_contigs[query_chr])
            else:
                self.query_sequence = self.query_contigs[query_chr]
            for chain in translocations:
                self.mash_fasta_and_chain_together(chain, False)


    def main(self, chromosome_name):
        import DDV

        names, query_chr, ref_chr = self.setup_for_reference_chromosome(chromosome_name)

        self.ref_sequence, contig_header = pluck_contig(ref_chr, names['ref'])  # only need the reference chromosome read, skip the others

        # Reference chain that establishes coordinate frame
        self.do_all_relevant_chains()

        names['query_gapped'], names['ref_gapped'] = self.write_gapped_fasta(names['query'], names['ref'])
        names['query_unique'], names['ref_unique'] = self.print_only_unique(names['query_gapped'], names['ref_gapped'])
        print("Finished creating gapped fasta files", names['ref'], names['query'])

        folder_name = self.output_folder_prefix + ref_chr
        source_path = '.\\bin\\Release\\output\\dnadata\\'
        if True:  #self.trial_run:  # these files are never used in the viz
            del names['ref']
            del names['query']
        self.move_fasta_source_to_destination(names, folder_name, source_path)
        if self.swap_columns:
            DDV.DDV_main(['DDV',
                          names['query_gapped'],
                          source_path, folder_name,
                          names['query_unique'],
                          names['ref_unique'],
                          names['ref_gapped']])
        else:
            DDV.DDV_main(['DDV',
                          names['ref_gapped'],
                          source_path, folder_name,
                          names['ref_unique'],
                          names['query_unique'],
                          names['query_gapped']])

    def setup_for_reference_chromosome(self, chromosome_name):
        if isinstance(chromosome_name, str):
            chromosome_name = (chromosome_name, chromosome_name)
        query_chr, ref_chr = chromosome_name
        names = {'ref': ref_chr + '_%s.fa' % first_word(self.ref_source),
                 'query': query_chr + '_%s.fa' % first_word(self.query_source)}  # for collecting all the files names in a modifiable way
        self.read_seq_to_memory(query_chr, ref_chr, self.query_source, self.ref_source, names['query'], names['ref'])

        self.relevant_chains = [chain for chain in self.chain_list if chain.tName == ref_chr]

        return names, query_chr, ref_chr


def do_chromosome(chr):
    parser = ChainParser(chain_name='hg38ToPanTro4.over.chain',
                         first_source='hg38_chr20.fa',  # 'HongKong\\hg38.fa',
                         second_source='panTro4.fa',  # 'panTro4_chr20.fa',
                         output_folder_prefix='panTro4_and_Hg38_trans1_',
                         trial_run=True,
                         swap_columns=True)
    # parser = ChainParser(chain_name='HongKong\\human_gorilla.bland.chain',
    #                      ref_source='HongKong\\susie3_agp.fasta',
    #                      query_source='HongKong\\hg38.fa',
    #                      output_folder_prefix='Susie3_and_Hg38_short3_',
    #                      trial_run=True,
    #                      swap_columns=True)
    parser.main(chr)


if __name__ == '__main__':
    chromosomes = [('chr20', 'chr20')]  # 'chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr22 chrX'.split()
    # TODO: handle Chr2A and Chr2B separately
    for chr in chromosomes:
        do_chromosome(chr)

    # import multiprocessing
    # workers = multiprocessing.Pool(6)  # number of simultaneous processes.  Watch your RAM usage
    # workers.map(do_chromosome, chromosomes)

