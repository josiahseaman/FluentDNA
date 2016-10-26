import os
from array import array
from datetime import datetime
import _io
from blist import blist

from ChainFiles import chain_file_to_list
from DDVUtils import just_the_name, chunks, pluck_contig, first_word, Batch, make_output_dir_with_suffix, ReverseComplement
from Span import AlignedSpans, Span


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
                 swap_columns=False, separate_translocations=False, squish_gaps=False,
                 show_translocations_only=False, aligned_only=False):
        self.ref_source = first_source  # example hg38ToPanTro4.chain  hg38 is the reference, PanTro4 is the query (has strand flips)
        self.query_source = second_source
        self.output_prefix = output_prefix
        self.output_folder = None
        self.query_contigs = dict()
        self.trial_run = trial_run
        self.swap_columns = swap_columns
        self.separate_translocations = separate_translocations
        self.show_translocations_only = show_translocations_only
        self.squish_gaps = squish_gaps
        self.aligned_only = aligned_only
        self.query_sequence = ''
        self.ref_sequence = ''
        self.query_seq_gapped = array('u', '')
        self.ref_seq_gapped = array('u', '')
        self.alignment = blist()  # optimized for inserts in the middle
        self.relevant_chains = []  # list of contigs that contribute to ref_chr in chain entries
        self.stored_rev_comps = {}
        self.output_fastas = []
        self.gapped = '_gapped'
        self.translocation_searched = 0
        self.translocation_deleted = 0

        self.chain_list = chain_file_to_list(chain_name)
        self.read_query_contigs(self.query_source)


    def read_query_contigs(self, input_file_path):
        print("Reading contigs... ", input_file_path)
        start_time = datetime.now()
        self.query_contigs = {}
        current_name = just_the_name(input_file_path)  # default to filename
        seq_collection = []

        # Pre-read generates an array of contigs with labels and sequences
        with open(input_file_path, 'r') as streamFASTAFile:
            for read in streamFASTAFile:
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
                    seq_collection.append(read.upper().rstrip())

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
        ref_pointer, query_pointer = self.process_chain_body(chain, ref_pointer, query_pointer, is_master_alignment)

        self.append_unaligned_end_in_master(chain, ref_pointer, query_pointer, is_master_alignment)


    def append_unaligned_end_in_master(self, chain, ref_pointer, query_pointer, is_master_alignment):
        if is_master_alignment and not self.trial_run:  # include unaligned ends
            ref_end = Span(ref_pointer, ref_pointer, chain.tName, chain.tStrand)
            query_end = Span(query_pointer, query_pointer, chain.qName, chain.qStrand)
            self.alignment.append(AlignedSpans(ref_end, query_end,
                                               len(self.query_sequence) - query_pointer,
                                               len(self.ref_sequence) - ref_pointer),)


    def alignment_chopping_index(self, new_alignment):
        """Return the index where to insert item x in list a, assuming a is sorted.

        The return value i is such that all e in a[:i] have e < x, and all e in
        a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
        insert just before the leftmost x already there.

        Optional args lo (default 0) and hi (default len(a)) bound the
        slice of a to be searched.
        """
        lo = 0
        hi = len(self.alignment)

        while lo < hi:
            mid = (lo + hi) // 2
            if self.alignment[mid] < new_alignment:
                lo = mid + 1
            else:
                hi = mid
        return lo


    def find_old_query_location(self, new_alignment, lo=0, hi=None, depth=0):
        """Binary search over the _mostly_ sorted list of query indices.
        Needs special casing to avoid looking at non-master chain entries."""
        if depth > 20:
            return False  # abandon hope
        if hi is None:
            self.translocation_searched += 1
        hi = hi or len(self.alignment)
        while lo < hi:
            mid = (lo + hi) // 2
            # Special handling for accidentally landing on a translocation (not sorted)
            if not self.alignment[mid].is_master_chain:  # recursively splits the binary search
                mid_left, mid_right = mid, mid
                # slide left
                while not self.alignment[mid_left].is_master_chain:
                    mid_left -= 1
                if not self.find_old_query_location(new_alignment, lo=lo, hi=mid_left, depth=depth + 1):
                    # slide right
                    while not self.alignment[mid_right].is_master_chain:
                        mid_right += 1
                    return self.find_old_query_location(new_alignment, lo=mid_right, hi=hi, depth=depth + 1)
                return True
            if self.alignment[mid].query_less_than(new_alignment):
                lo = mid + 1
            else:
                hi = mid

        lo = max(lo - 1, 0)
        final_possible = self.alignment[lo].query_unique_span()
        while new_alignment.query.begin not in final_possible and new_alignment.query.end >= final_possible.end:# and lo < hi:
            lo += 1
            final_possible = self.alignment[lo].query_unique_span()
        if new_alignment.query.begin not in final_possible or new_alignment.query.end not in final_possible:
            return False  # the old query copy is being used in more than one alignment pair
        #     # raise ValueError("%s   %s   %s" % tuple(str(self.alignment[x].query_unique_span()) for x in [lo-1, lo, lo+1]))
        self.alignment[lo].remove_old_query_copy(new_alignment)
        self.translocation_deleted += 1

        print(int(self.translocation_deleted / self.translocation_searched * 100), "%", self.translocation_searched, self.translocation_deleted )
        return True


    def process_chain_body(self, chain, ref_pointer, query_pointer, is_master_alignment):
        assert len(chain.entries), "Chain has no data"
        for entry_index, entry in enumerate(chain.entries):  # ChainEntry
            size, gap_query, gap_ref = entry.size, entry.gap_query, entry.gap_ref
            first_in_chain = entry_index == 0
            if not size:
                continue  # entries with 0 size don't count
            # Debugging code
            if is_master_alignment and self.trial_run and len(self.alignment) > 1000:  # 9500000  is_master_alignment and
                break
            if not is_master_alignment:  # skip the unaligned middle of translocation chains
                if self.separate_translocations and max(gap_query, gap_ref) > 2600:
                    first_in_chain = True
                    gap_ref, gap_query = 0, 0  # This is caused by overzealous netting
                elif not self.separate_translocations:
                    gap_ref, gap_query = 0, 0  # placed translocations don't need gaps

            aligned_query = Span(query_pointer, query_pointer + size, chain.qName, chain.qStrand)
            aligned_ref = Span(ref_pointer, ref_pointer + size, chain.tName, chain.tStrand)
            new_alignment = AlignedSpans(aligned_ref, aligned_query, gap_ref, gap_query, is_master_alignment, first_in_chain)

            if is_master_alignment or self.separate_translocations:
                self.alignment.append(new_alignment)
            else:
                # Remove old unaligned query location
                scrutiny_index = self.find_old_query_location(new_alignment)  # Binary search using query


                # Add new_alignment at ref location
                scrutiny_index = max(self.alignment_chopping_index(new_alignment) - 1, 0)  # Binary search
                old = self.alignment.pop(scrutiny_index)
                first, second = old.align_ref_unique(new_alignment)
                self.alignment.insert(scrutiny_index, first)
                self.alignment.insert(scrutiny_index + 1, second)

            query_pointer += size + entry.gap_ref  # alignable and unalignable block concatenated together
            ref_pointer += size + entry.gap_query  # two blocks of sequence separated by gap

            # TODO handle interlacing
        return ref_pointer, query_pointer


    def create_fasta_from_composite_alignment(self):
        """self.alignment is a data structure representing the composites of all the relevant
        chain data.  This method turns that data construct into a gapped FASTA file by reading the original
        FASTA files."""
        previous_chr = None
        for pair in self.alignment:
            if previous_chr != (pair.query.contig_name, pair.query.strand):
                if not self.switch_sequences(pair.query.contig_name, pair.query.strand):  # pair.ref.contig_name could be None
                    raise FileNotFoundError("Could not switch to " + pair.query.contig_name)  # We'll have to skip this alignment because we don't have the FASTA for it
            previous_chr = (pair.query.contig_name, pair.query.strand)

            query_snippet = pair.query.sample(self.query_sequence)
            if not self.aligned_only:
                query_snippet += pair.query_unique_span().sample(self.query_sequence)
                query_snippet += 'X' * pair.ref_tail_size  # whenever there is no alignable sequence, it's filled with X's

            ref_snippet = pair.ref.sample(self.ref_sequence)
            if not self.aligned_only:  # Algined_only simply skips over the unaligned tails
                ref_snippet += 'X' * pair.query_tail_size  # Ref 'X' gap is in the middle, query is at the end, to alternate
                ref_snippet += pair.ref_unique_span().sample(self.ref_sequence)
            if self.show_translocations_only and pair.is_master_chain:  # main chain
                ref_snippet = 'X' * len(ref_snippet)
                query_snippet = 'X' * len(query_snippet)
            elif self.separate_translocations and pair.is_first_entry and not pair.is_master_chain:
                self.add_translocation_header(pair)
            self.query_seq_gapped.extend(query_snippet)
            self.ref_seq_gapped.extend(ref_snippet)


    def switch_sequences(self, query_name, query_strand):
        # TODO: self.ref_sequence = self.ref_contigs[ref_name]
        if query_name in self.query_contigs:
            if query_strand == '-':  # need to load rev_comp
                if query_name not in self.stored_rev_comps:
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
            self.alignment.append(AlignedSpans(Span(0, 0, chain.tName, chain.tStrand),
                                               Span(0, 0, chain.qName, chain.qStrand),
                                               query_pointer,
                                               ref_pointer))
        return ref_pointer, query_pointer


    def pad_next_line(self):
        self.ref_seq_gapped.extend('X' * len(self.ref_seq_gapped) % 100)
        self.query_seq_gapped.extend('X' * len(self.query_seq_gapped) % 100)


    def add_translocation_header(self, alignment):
        """ :param alignment: AlignedSpan
        """
        self.ref_seq_gapped.extend('\n>%s_%s_%i\n' % (alignment.ref.contig_name, alignment.ref.strand, alignment.ref.begin))  # visual separators
        self.query_seq_gapped.extend('\n>%s_%s_%i\n' % (alignment.query.contig_name, alignment.query.strand, alignment.query.begin))

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
            self.mash_fasta_and_chain_together(chain, is_master_alignment)  # first chain is the master alignment
            previous = chain
        self.create_fasta_from_composite_alignment()


    def setup_for_reference_chromosome(self, ref_chr):
        ending = ref_chr + '__squished' * self.squish_gaps + \
            '__separate_translocations' * self.separate_translocations + \
            '__translocations' * self.show_translocations_only + \
            '__aligned_only' * self.aligned_only
        self.output_folder = make_output_dir_with_suffix(self.output_prefix, ending)
        names = {'ref': ref_chr + '_%s.fa' % first_word(self.ref_source),
                 'query': '%s_to_%s_%s.fa' % (first_word(self.query_source), first_word(self.ref_source), ref_chr)
                 }  # for collecting all the files names in a modifiable way
        # This assumes the chains have been sorted by score, so the highest score is the matching query_chr
        self.relevant_chains = [chain for chain in self.chain_list if chain.tName == ref_chr]

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

