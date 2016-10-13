import os
#from typing import List

from array import array
from collections import defaultdict
from datetime import datetime

from ChainFiles import chain_file_to_list
from DDVUtils import just_the_name, chunks, pluck_contig, first_word, Batch, make_output_dir_with_suffix, ReverseComplement


class ChainParser:
    def __init__(self, chain_name, second_source, first_source, output_prefix, trial_run=False,
                 swap_columns=False, include_translocations=True, squish_gaps=False):
        self.width_remaining = defaultdict(lambda: 70)
        self.ref_source = first_source  # example hg38ToPanTro4.chain  hg38 is the reference, PanTro4 is the query (has strand flips)
        self.query_source = second_source
        self.output_prefix = output_prefix
        self.output_folder = None
        self.query_contigs = dict()
        self.trial_run = trial_run
        self.swap_columns = swap_columns
        self.include_translocations = include_translocations
        self.squish_gaps = squish_gaps
        self.query_sequence = ''
        self.ref_sequence = ''
        self.query_seq_gapped = array('u', '')
        self.ref_seq_gapped = array('u', '')
        self.relevant_chains = []  # list of contigs that contribute to ref_chr in chain entries
        self.stored_rev_comps = {}
        self.output_fastas = []
        self.gapped = '_gapped'

        self.chain_list = chain_file_to_list(chain_name)
        if self.include_translocations:
            self.read_contigs(self.query_source)


    def read_contigs(self, input_file_path):
        print("Reading contigs... ", input_file_path)
        start_time = datetime.now()
        self.query_contigs = {}
        current_name = just_the_name(input_file_path)  # default to filename
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
        print("Read %i FASTA Contigs in:" % len(self.query_contigs), datetime.now() - start_time)


    def read_query_seq_to_memory(self, query_chr, query_source):
        if not self.include_translocations:
            self.query_contigs[query_chr] = pluck_contig(query_chr, query_source)  # skip others
        else:
            pass  # already covered in __init__ query_contigs and do_all_relevant_chains setup


    def write_fasta_lines(self, filestream, seq):
        if isinstance(filestream, str):  # I'm actually given a file name and have to open it myself
            file_name = os.path.join(self.output_folder, filestream)
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
        assert len(chain.entries), "Chain has no data"
        for entry_index, entry in enumerate(chain.entries):  # ChainEntry
            size, gap_query, gap_ref = entry.size, entry.gap_query, entry.gap_ref

            # Debugging code
            if is_master_alignment and self.trial_run and len(self.ref_seq_gapped) > 1000000:  # 9500000  is_master_alignment and
                break

            space_saved = 0
            if not is_master_alignment and max(gap_query, gap_ref) > 2600:  # 26 lines is the height of a label
                # Don't show intervening sequence
                # #14 Skipping display of large gaps formed by bad netting.
                # TODO insert header
                gap_query, gap_ref = 0, 0
                # Pointer += uses unmodified entry.gap_ref, so skips the full sequence
            elif self.squish_gaps:
                space_saved = min(gap_query, gap_ref)

            query_seq_absolute = self.query_sequence[query_pointer: query_pointer + size + gap_ref]
            query_snippet = query_seq_absolute + 'X' * (gap_query - space_saved)
            query_pointer += size + entry.gap_ref  # alignable and unalignable block concatenated together
            self.query_seq_gapped.extend(query_snippet)

            ref_snippet = self.ref_sequence[ref_pointer: ref_pointer + size] + 'X' * (gap_ref - space_saved)
            ref_snippet += self.ref_sequence[ref_pointer + size: ref_pointer + size + gap_query]
            ref_pointer += size + entry.gap_query  # two blocks of sequence separated by gap
            self.ref_seq_gapped.extend(ref_snippet)

            if len(ref_snippet) != len(query_snippet):
                if len(ref_snippet) < len(query_snippet):
                    faulty_source, faulty_chr = self.ref_source, chain.tName
                else:
                    faulty_source, faulty_chr = self.query_source, chain.qName
                remaining_entries = len(chain.entries) - entry_index
                raise EOFError("Reached the end of %s before expected point.\n"
                               "Chain file still has %i entries on %s to process." % (faulty_source, remaining_entries, faulty_chr))

        if is_master_alignment and not self.trial_run:  # last one: print out all remaining sequence
            gap_query, gap_ref = len(self.ref_sequence) - ref_pointer, len(self.query_sequence) - query_pointer
            self.query_seq_gapped.extend(self.query_sequence[query_pointer:] + 'X' * gap_query)
            self.ref_seq_gapped.extend('X' * gap_ref + self.ref_sequence[ref_pointer:])


    def setup_chain_start(self, chain, is_master_alignment):
        # convert to int except for chr names and strands
        ref_pointer, query_pointer = chain.tStart, chain.qStart
        if chain.tEnd - chain.tStart > 100 * 1000:
            print('>>>>', chain)
        if is_master_alignment and not self.trial_run:  # include the unaligned beginning of the sequence
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
        query_gap_name = os.path.join(self.output_folder, os.path.splitext(query)[0] + self.gapped + '.fa')
        ref_gap_name = os.path.join(self.output_folder, os.path.splitext(reference)[0] + self.gapped + '.fa')
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
        query_unique_name = os.path.join(self.output_folder, query_gapped_name.replace('%s' % self.gapped, '_unique'))
        ref_unique_name = os.path.join(self.output_folder, ref_gapped_name.replace(self.gapped, '_unique'))
        with open(query_unique_name, 'w') as query_filestream:
            self.write_fasta_lines(query_filestream, ''.join(query_uniq_array))
        with open(ref_unique_name, 'w') as ref_filestream:
            self.write_fasta_lines(ref_filestream, ''.join(que_uniq_array))

        return query_unique_name, ref_unique_name


    def do_all_relevant_chains(self, relevant_chains=None):
        """In panTro4ToHg38.over.chain there are ZERO chains that have a negative strand on the reference 'tStrand'.
        I think it's a rule that you always flip the query strand instead."""
        if relevant_chains is None:
            relevant_chains = self.relevant_chains
        previous = None
        for chain in relevant_chains:
            is_master_alignment = previous is None
            if self.include_translocations or is_master_alignment:  # otherwise we'll only do the master alignment
                if chain.qName in self.query_contigs:
                    if chain.qStrand == '-':  # need to load rev_comp
                        if chain.qName not in self.stored_rev_comps:
                            if len(self.query_contigs[chain.qName]) > 1000000:
                                print("Reversing", chain.qName, len(self.query_contigs[chain.qName]) // 1000000, 'Mbp')
                            self.stored_rev_comps[chain.qName] = ReverseComplement(self.query_contigs[chain.qName])  # caching for performance
                            # TODO: replace rev_comp with a lazy generator that mimics having the whole string using [x] and [y:x]
                        self.query_sequence = self.stored_rev_comps[chain.qName]
                    else:
                        self.query_sequence = self.query_contigs[chain.qName]
                    self.mash_fasta_and_chain_together(chain, is_master_alignment)  # first chain is the master alignment
                    previous = chain
                else:
                    print("No fasta source for", chain.qName)


    def setup_for_reference_chromosome(self, chromosome_name):
        query_chr, ref_chr = chromosome_name, chromosome_name
        ending = ref_chr + '__squished' * self.squish_gaps + '__no_translocations' * (not self.include_translocations)
        self.output_folder = make_output_dir_with_suffix(self.output_prefix, ending)
        names = {'ref': ref_chr + '_%s.fa' % first_word(self.ref_source),
                 'query': '%s_to_%s_%s.fa' % (first_word(self.query_source), first_word(self.ref_source), query_chr)
                 }  # for collecting all the files names in a modifiable way
        self.read_query_seq_to_memory(query_chr, self.query_source)

        self.relevant_chains = [chain for chain in self.chain_list if chain.tName == ref_chr]

        return names, query_chr, ref_chr


    def _parse_chromosome_in_chain(self, chromosome_name) -> Batch:
        names, query_chr, ref_chr = self.setup_for_reference_chromosome(chromosome_name)

        self.ref_sequence = pluck_contig(ref_chr, self.ref_source)  # only need the reference chromosome read, skip the others

        self.do_all_relevant_chains()

        names['query_gapped'], names['ref_gapped'] = self.write_gapped_fasta(names['query'], names['ref'])
        names['query_unique'], names['ref_unique'] = self.print_only_unique(names['query_gapped'], names['ref_gapped'])
        # NOTE: Order of these appends DOES matter!
        self.output_fastas.append(names['ref_gapped'])
        self.output_fastas.append(names['ref_unique'])
        self.output_fastas.append(names['query_unique'])
        self.output_fastas.append(names['query_gapped'])
        if self.swap_columns:
            self.output_fastas = self.output_fastas.reverse()
        print("Finished creating gapped fasta files", names['ref'], names['query'])

        if True:  #self.trial_run:  # these files are never used in the viz
            del names['ref']
            del names['query']
        batch = Batch(chromosome_name, self.output_fastas, self.output_folder)
        self.output_folder = None  # clear the previous value
        return batch


    def parse_chain(self, chromosomes) -> list:
        assert isinstance(chromosomes, list), "'Chromosomes' must be a list! A single element list is okay."

        batches = []
        for chromosome in chromosomes:
            batches.append(self._parse_chromosome_in_chain(chromosome))
        return batches
        # workers = multiprocessing.Pool(6)  # number of simultaneous processes. Watch your RAM usage
        # workers.map(self._parse_chromosome_in_chain, chromosomes)
