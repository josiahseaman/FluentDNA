import os
from array import array
from datetime import datetime
import _io
from blist import blist

from ChainFiles import chain_file_to_list
from DDVUtils import just_the_name, chunks, pluck_contig, first_word, Batch, make_output_dir_with_suffix, ReverseComplement
from Span import AlignedPair, Span


def scan_past_header(seq, index):
    """Moves the pointer past any headers or comments and places index on the next valid sequence character.
    Doesn't increment at all if it is called in a place that is not the start of a header."""
    if seq[index] == '\n':  # skip newline marking the end of a contig
        index += 1
    scanning_past_header = seq[index] in '>;'  # ; is for comments
    while scanning_past_header:
        scanning_past_header = seq[index] != '\n'
        index += 1  # skip this character

    return index


class ChainParser:
    def __init__(self, chain_name, second_source, first_source, output_prefix, trial_run=False,
                 swap_columns=False, include_translocations=True, squish_gaps=False):
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
        self.alignment = blist()  # optimized for inserts in the middle
        self.relevant_chains = []  # list of contigs that contribute to ref_chr in chain entries
        self.stored_rev_comps = {}
        self.output_fastas = []
        self.gapped = '_gapped'

        self.chain_list = chain_file_to_list(chain_name)
        if self.include_translocations:
            self.read_query_contigs(self.query_source)


    def read_query_contigs(self, input_file_path):
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



    def _write_fasta_lines(self, filestream, seq):
        assert isinstance(filestream, _io.TextIOWrapper)  # I'm actually given a file name and have to open it myself
        contigs = seq.split('\n')
        index = 0
        while index < len(contigs):
            if contigs[index].startswith('>'):
                header, contents = contigs[index], contigs[index + 1]
                index += 2
            else:
                header, contents = None, contigs[index]
                index += 1
            self.__do_write(filestream, contents, header)


    def __do_write(self, filestream, seq, header=None):
        """Specialized function for writing sets of headers and sequence in FASTA.
        It chunks the file up into 70 character lines, but leaves headers alone"""
        if header is not None:
            filestream.write(header + '\n')  # double check newlines
        try:
            for line in chunks(seq, 70):
                filestream.write(line + '\n')
        except Exception as e:
            print(e)


    def mash_fasta_and_chain_together(self, chain, is_master_alignment=False):
        ref_pointer, query_pointer = self.setup_chain_start(chain, is_master_alignment)

        # if not is_master_alignment:
        #     self.do_translocation_housework(chain, query_pointer, ref_pointer)
        self.process_chain_body(chain, query_pointer, ref_pointer, is_master_alignment)


    def process_chain_body(self, chain, query_pointer, ref_pointer, is_master_alignment):
        assert len(chain.entries), "Chain has no data"
        for entry_index, entry in enumerate(chain.entries):  # ChainEntry
            size, gap_query, gap_ref = entry.size, entry.gap_query, entry.gap_ref

            # Debugging code
            if is_master_alignment and self.trial_run and len(self.ref_seq_gapped) > 1000000:  # 9500000  is_master_alignment and
                break

            aligned_query = Span(query_pointer, query_pointer + size, chain.qName, chain.qStrand)
            aligned_ref = Span(ref_pointer, ref_pointer + size, chain.tName, chain.tStrand)
            self.alignment.append(AlignedPair(aligned_ref, aligned_query))
            self.alignment.append(AlignedPair(None, Span(query_pointer + size, query_pointer + size + gap_ref, chain.qName, chain.qStrand)))
            self.alignment.append(AlignedPair(Span(ref_pointer + size, ref_pointer + size + gap_query, chain.tName, chain.tStrand), None))
            query_pointer += size + entry.gap_ref  # alignable and unalignable block concatenated together
            ref_pointer += size + entry.gap_query  # two blocks of sequence separated by gap

            # TODO handle interlacing


    def create_fasta_from_composite_alignment(self):
        previous_chr = None
        for pair in self.alignment:
            if previous_chr != (pair.query.contig_name, pair.query.strand):
                if not self.switch_sequences(pair.ref.contig_name, pair.query.contig_name, pair.query.strand):
                    continue  # We'll have to skip this alignment because we don't have the FASTA for it
                previous_chr = (pair.query.contig_name, pair.query.strand)
            if pair.query is None:  # whenever there is no alignable sequence, it's filled with N's
                query_snippet = 'X' * pair.ref.size()
            else:
                query_snippet = self.query_sequence[pair.query.begin: pair.query.end]
            self.query_seq_gapped.extend(query_snippet)
            if pair.ref is None:
                ref_snippet = 'X' * pair.query.size()
            else:
                ref_snippet = self.ref_sequence[pair.ref.begin: pair.ref.end]
            self.ref_seq_gapped.extend(ref_snippet)


    def switch_sequences(self, ref_name, query_name, query_strand):
        # TODO: self.ref_sequence = self.ref_contigs[ref_name]
        if query_name in self.query_contigs:
            if query_strand == '-':  # need to load rev_comp
                if query_name not in self.stored_rev_comps:
                    if len(self.query_contigs[query_name]) > 1000000:
                        print("Reversing", query_name, len(self.query_contigs[query_name]) // 1000000, 'Mbp')
                    self.stored_rev_comps[query_name] = ReverseComplement(self.query_contigs[query_name])  # caching for performance
                self.query_sequence = self.stored_rev_comps[query_name]
            else:
                self.query_sequence = self.query_contigs[query_name]
        else:
            print("ERROR: No fasta source for", query_name)
            return False
        return True


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


    def pad_next_line(self):
        self.ref_seq_gapped.extend('X' * len(self.ref_seq_gapped) % 100)
        self.query_seq_gapped.extend('X' * len(self.query_seq_gapped) % 100)


    def do_translocation_housework(self, chain, query_pointer, ref_pointer):
        self.ref_seq_gapped.extend('\n>%s_%s_%i\n' % (chain.tName, chain.tStrand, ref_pointer))  # visual separators
        self.query_seq_gapped.extend('\n>%s_%s_%i\n' % (chain.qName, chain.qStrand, query_pointer))

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
        self.write_complete_fasta(query_gap_name, self.query_seq_gapped)
        self.write_complete_fasta(ref_gap_name, self.ref_seq_gapped)
        return query_gap_name, ref_gap_name


    def write_complete_fasta(self, file_path, seq_content_array):
        """This function ensures that all FASTA files start with a >header\n line"""
        with open(file_path, 'w') as filestream:
            if seq_content_array[0] != '>':  # start with a header
                temp_content = seq_content_array
                header = '>%s\n' % just_the_name(file_path)
                if isinstance(temp_content, list):
                    seq_content_array = [header]
                else:
                    seq_content_array = array('u', header)
                seq_content_array.extend(temp_content)
            self._write_fasta_lines(filestream, ''.join(seq_content_array))


    def print_only_unique(self, query_gapped_name, ref_gapped_name):
        query_uniq_array = array('u', self.query_seq_gapped)
        ref_uniq_array = array('u', self.ref_seq_gapped)
        print("Done allocating unique array")
        r, q = 0, 0  # indices
        while q < len(self.query_seq_gapped) and r < len(self.ref_seq_gapped):
            q = scan_past_header(self.query_seq_gapped, q)  # query_uniq_array is already initialized to contain header characters
            r = scan_past_header(self.ref_seq_gapped, r)

            # only overlapping section
            if self.query_seq_gapped[q] == self.ref_seq_gapped[r]:
                query_uniq_array[q] = 'X'
                ref_uniq_array[r] = 'X'
            else:  # Not equal
                if self.query_seq_gapped[q] == 'N':
                    query_uniq_array[q] = 'X'
                if self.ref_seq_gapped[r] == 'N':
                    ref_uniq_array[r] = 'X'
            q += 1
            r += 1

        # Just to be thorough: prints aligned section (shortest_sequence) plus any dangling end sequence
        query_unique_name = os.path.join(self.output_folder, query_gapped_name.replace(self.gapped, '_unique'))
        ref_unique_name = os.path.join(self.output_folder, ref_gapped_name.replace(self.gapped, '_unique'))
        self.write_complete_fasta(query_unique_name, query_uniq_array)
        self.write_complete_fasta(ref_unique_name, ref_uniq_array)

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
                self.mash_fasta_and_chain_together(chain, is_master_alignment)  # first chain is the master alignment
                previous = chain
        self.create_fasta_from_composite_alignment()


    def setup_for_reference_chromosome(self, ref_chr):
        ending = ref_chr + '__squished' * self.squish_gaps + '__no_translocations' * (not self.include_translocations)
        self.output_folder = make_output_dir_with_suffix(self.output_prefix, ending)
        names = {'ref': ref_chr + '_%s.fa' % first_word(self.ref_source),
                 'query': '%s_to_%s_%s.fa' % (first_word(self.query_source), first_word(self.ref_source), ref_chr)
                 }  # for collecting all the files names in a modifiable way
        # This assumes the chains have been sorted by score, so the highest score is the matching query_chr
        self.relevant_chains = [chain for chain in self.chain_list if chain.tName == ref_chr]
        if self.include_translocations and self.query_source:
            query_chr = self.relevant_chains[0].qName
            self.query_contigs[query_chr] = pluck_contig(query_chr, self.query_source)
        else:  # read in all the contigs
            pass  # already covered in __init__ query_contigs and do_all_relevant_chains setup
        return names, ref_chr


    def _parse_chromosome_in_chain(self, chromosome_name) -> Batch:
        names, ref_chr = self.setup_for_reference_chromosome(chromosome_name)

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

