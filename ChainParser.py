import os
from datetime import datetime

import blist
from array import array

from ChainFiles import chain_file_to_list, match
from DDVUtils import make_output_dir_with_suffix
from DefaultOrderedDict import DefaultOrderedDict
from Span import AlignedSpans, Span, alignment_chopping_index
from DNASkittleUtils.CommandLineUtils import just_the_name
from DNASkittleUtils.DDVUtils import first_word, Batch,  ReverseComplement, BlankIterator
from DNASkittleUtils.Contigs import pluck_contig, write_complete_fasta

def scan_past_header(seq, index, take_shortcuts=False, skip_newline=True):
    """Moves the pointer past any headers or comments and places index on the next valid sequence character.
    Doesn't increment at all if it is called in a place that is not the start of a header.
    NOTE: this will crash with IndexError if there's a header at the end of the file with no body."""
    if take_shortcuts:
        return index
    if seq[index] == '\n':  # skip newline marking the end of a contig
        index += 1
    if seq[index] not in '>;':  # ; is for comments
        return index
    while seq[index] != '\n':
        index += 1  # skip this character
    if seq[index] == '\n' and skip_newline:  # skip newline marking the end of a header
        index += 1

    return index


class ChainParser:
    def __init__(self, chain_name, first_source, second_source, output_prefix,
                 trial_run=False, separate_translocations=False, squish_gaps=False,
                 show_translocations_only=False, aligned_only=False, no_titles=False):
        self.ref_source = first_source  # example hg38ToPanTro4.chain  hg38 is the reference, PanTro4 is the query (has strand flips)
        self.query_source = second_source
        self.output_prefix = output_prefix
        self.output_folder = None
        self.query_contigs = dict()
        self.trial_run = trial_run
        self.separate_translocations = separate_translocations
        self.no_titles = no_titles
        self.show_translocations_only = show_translocations_only
        self.squish_gaps = squish_gaps
        self.aligned_only = aligned_only
        self.query_sequence = ''
        self.ref_sequence = ''
        self.query_seq_gapped = array('u', '')
        self.ref_seq_gapped = array('u', '')
        self.output_fastas = []
        self.alignment = blist.blist()  # optimized for inserts in the middle
        self.stored_rev_comps = {}
        self.gapped = '_gapped'
        self.stats = DefaultOrderedDict(lambda: 0)

        self.read_query_contigs(self.query_source)
        self.chain_list = chain_file_to_list(chain_name)


    def write_stats_file(self):
        with open(os.path.join(self.output_folder, 'stats.txt'), 'w+') as stats:
            stats.write('\n===== Alignment Stats ======\n')
            for key in self.stats:
                stats.write('%s\t%i\n' % (key, self.stats[key]))
            stats.write('\n====================================\n')


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


    def mash_fasta_and_chain_together(self, chain, is_master_alignment=False):
        ref_pointer, query_pointer = self.setup_chain_start(chain, is_master_alignment)
        ref_pointer, query_pointer = self.process_chain_body(chain, ref_pointer, query_pointer, is_master_alignment)
        self.append_unaligned_end_in_master(chain, ref_pointer, query_pointer, is_master_alignment)


    def append_unaligned_end_in_master(self, chain, ref_pointer, query_pointer, is_master_alignment):
        if is_master_alignment:
            if not self.query_sequence:
                self.switch_sequences(chain.qName, chain.qStrand)
            if not self.trial_run and self.query_sequence and self.ref_sequence:  # include unaligned ends
                ref_end = Span(ref_pointer, ref_pointer, chain.tName, chain.tStrand)
                query_end = Span(query_pointer, query_pointer, chain.qName, chain.qStrand)
                self.alignment.append(AlignedSpans(ref_end, query_end,
                                                   len(self.query_sequence) - query_pointer,
                                                   len(self.ref_sequence) - ref_pointer),)

    def find_old_query_location(self, new_alignment, lo=0, hi=None, depth=0):
        """Binary search over the _mostly_ sorted list of query indices.
        Needs special casing to avoid looking at non-master chain entries."""
        if depth > 30:
            return False  # abandon hope
        if hi is None:
            self.stats['translocation_searched'] += 1
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
            # Actual Binary search is here:
            if new_alignment.query_less_than(self.alignment[mid]):
                hi = mid
            else:
                lo = mid + 1

        # Binary search brings us to the right most position
        lo = max(lo - 1, 0)
        final_possible = self.alignment[lo].query_unique_span()
        # while new_alignment.query.begin not in final_possible and new_alignment.query.end >= final_possible.end:# and lo < hi:
        #     lo += 1
        #     final_possible = self.alignment[lo].query_unique_span()
        if new_alignment.query.begin not in final_possible or new_alignment.query.end not in final_possible:
            return False  # the old query copy is being used in more than one alignment pair
        #     # raise ValueError("%s   %s   %s" % tuple(str(self.alignment[x].query_unique_span()) for x in [lo-1, lo, lo+1]))
        replacements = self.alignment.pop(lo).remove_old_query_copy(new_alignment)
        [self.alignment.insert(lo, x) for x in reversed(replacements)]  # insert them into alignment in place, in order
        self.stats['translocation_deleted'] += 1

        print(int(self.stats['translocation_deleted'] / self.stats['translocation_searched'] * 100),
              "%", self.stats['translocation_searched'], self.stats['translocation_deleted'])
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
                # self.find_old_query_location(new_alignment)  # Binary search using query

                # Add new_alignment at ref location
                scrutiny_index = max(alignment_chopping_index(self.alignment, new_alignment) - 1, 0)  # Binary search
                try:
                    old = self.alignment.pop(scrutiny_index)
                    first, second = old.align_ref_unique(new_alignment)
                    self.alignment.insert(scrutiny_index, first)
                    self.alignment.insert(scrutiny_index + 1, second)
                except IndexError as e:
                    print(e)

            if entry.gap_query > 0: self.stats['Query Number of Gaps (all)'] += 1
            if entry.gap_query > 10: self.stats['Query Gaps larger than 10bp'] += 1
            if entry.gap_query > 100: self.stats['Query Gaps larger than 100bp'] += 1
            if entry.gap_query > 1000: self.stats['Query Gaps larger than 1000bp'] += 1

            if entry.gap_ref > 0: self.stats['Ref Number of Gaps (all)'] += 1
            if entry.gap_ref > 10: self.stats['Ref Gaps larger than 10bp'] += 1
            if entry.gap_ref > 100: self.stats['Ref Gaps larger than 100bp'] += 1
            if entry.gap_ref > 1000: self.stats['Ref Gaps larger than 1000bp'] += 1
            query_pointer += size + entry.gap_ref  # alignable and unalignable block concatenated together
            ref_pointer += size + entry.gap_query  # two blocks of sequence separated by gap

            # TODO handle interlacing
        return ref_pointer, query_pointer


    def create_fasta_from_composite_alignment(self, previous_chr=None, alignment=None):
        """self.alignment is a data structure representing the composites of all the relevant
        chain data.  This method turns that data construct into a gapped FASTA file by reading the original
        FASTA files."""
        if alignment is None:
            alignment = self.alignment

        for pair in alignment:
            if previous_chr != (pair.query.contig_name, pair.query.strand):
                if not self.switch_sequences(pair.query.contig_name, pair.query.strand):  # pair.ref.contig_name could be None
                    continue  # skip this pair since it can't be displayed
            previous_chr = (pair.query.contig_name, pair.query.strand)

            query_snippet = pair.query.sample(self.query_sequence)
            if not self.aligned_only:
                query_snippet += pair.query_unique_span().sample(self.query_sequence)
                query_snippet += 'X' * pair.ref_tail_size  # whenever there is no alignable sequence, it's filled with X's

            ref_snippet = pair.ref.sample(self.ref_sequence)
            if not self.aligned_only:  # Aligned_only simply skips over the unaligned tails
                ref_snippet += 'X' * pair.query_tail_size  # Ref 'X' gap is in the middle, query is at the end, to alternate
                ref_snippet += pair.ref_unique_span().sample(self.ref_sequence)

            if pair.is_hidden or self.show_translocations_only and pair.is_master_chain:  # main chain
                ref_snippet = 'X' * len(ref_snippet)
                query_snippet = 'X' * len(query_snippet)
            elif self.separate_translocations and pair.is_first_entry and not pair.is_master_chain:
                self.add_translocation_header(pair)
            if len(query_snippet) != len(ref_snippet):  # make absolute sure we don't step out of phase with bad lengths
                difference = len(ref_snippet) - len(query_snippet)
                ref_snippet += 'X' * max(0, -difference)
                query_snippet += 'X' * max(0, difference)

            assert len(query_snippet) == len(ref_snippet), "Mismatched lengths: " + '\n'.join([ref_snippet, query_snippet])
            self.query_seq_gapped.extend(query_snippet)
            self.ref_seq_gapped.extend(ref_snippet)

            # else:  # if 'random' not in str(pair) and 'unknown' not in str(pair):  # Only notify me when a
            #     print("Processed bad alignment pair:\n", pair, "\nRef:", ref_snippet, "\nQry:", query_snippet)

        print("Done gapping sequence")


    def switch_sequences(self, query_name, query_strand):
        # TODO: self.ref_sequence = self.ref_contigs[ref_name]
        if query_name in self.query_contigs:  # TODO: capitalization!!!!
            if query_strand == '-':  # need to load rev_comp
                if query_name not in self.stored_rev_comps:
                    self.stored_rev_comps[query_name] = self.rev_comp_contig(query_name)  # caching for performance
                self.query_sequence = self.stored_rev_comps[query_name]
            else:
                self.query_sequence = self.query_contigs[query_name]
        else:
            return self.missing_query_sequence(query_name)
        return True


    def missing_query_sequence(self, query_name):
        print("ERROR: No fasta source for", query_name)
        self.query_sequence = BlankIterator('N')
        return False


    def rev_comp_contig(self, query_name):
        return ReverseComplement(self.query_contigs[query_name])


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
        column_width = 100
        characters_remaining = column_width - (len(self.ref_seq_gapped) % column_width)
        self.ref_seq_gapped.extend('X' * characters_remaining)
        self.query_seq_gapped.extend('X' * characters_remaining)


    def add_translocation_header(self, alignment):
        """ :param alignment: AlignedSpan
        """
        if self.no_titles:
            self.pad_next_line()
        else:
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


    def write_gapped_fasta(self, reference, query, prepend_output_folder=True):
        query_gap_name = os.path.splitext(query)[0] + self.gapped + '.fa'
        ref_gap_name = os.path.splitext(reference)[0] + self.gapped + '.fa'
        if prepend_output_folder:
            query_gap_name = os.path.join(self.output_folder, query_gap_name)
            ref_gap_name = os.path.join(self.output_folder, ref_gap_name)
        write_complete_fasta(query_gap_name, self.query_seq_gapped)
        write_complete_fasta(ref_gap_name, self.ref_seq_gapped)
        print("Finished creating gapped fasta files", ref_gap_name, query_gap_name)
        return ref_gap_name, query_gap_name


    def print_only_unique(self, query_gapped_name, ref_gapped_name):
        query_uniq_array, ref_uniq_array = self.compute_unique_sequence()

        query_unique_name = query_gapped_name.replace(self.gapped, '_unique')
        if not query_unique_name.startswith(self.output_folder):  # TODO: this shouldn't happen but it does
            query_unique_name = os.path.join(self.output_folder, query_unique_name)
        ref_unique_name = ref_gapped_name.replace(self.gapped, '_unique')
        if not ref_unique_name.startswith(self.output_folder):
            ref_unique_name = os.path.join(self.output_folder, ref_unique_name)
        write_complete_fasta(query_unique_name, query_uniq_array)
        write_complete_fasta(ref_unique_name, ref_uniq_array)

        return ref_unique_name, query_unique_name

    def compute_unique_sequence(self):
        query_uniq_array = array('u', self.query_seq_gapped)
        ref_uniq_array = array('u', self.ref_seq_gapped)
        print("Done allocating unique array")
        # query_uniq_array is already initialized to contain header characters
        q = scan_past_header(self.query_seq_gapped, 0)
        r = scan_past_header(self.ref_seq_gapped, 0)
        while q < len(self.query_seq_gapped) and r < len(self.ref_seq_gapped):
            q = scan_past_header(self.query_seq_gapped, q)  # query_uniq_array is already initialized to contain header characters
            r = scan_past_header(self.ref_seq_gapped, r)

            # only overlapping section
            q_letter = self.query_seq_gapped[q]
            r_letter = self.ref_seq_gapped[r]
            if q_letter == r_letter:
                query_uniq_array[q] = 'X'
                ref_uniq_array[r] = 'X'
                self.stats['Shared seq bp'] += 1
            else:  # Not equal
                if q_letter == 'N':
                    query_uniq_array[q] = 'X'
                    self.stats['Query N to ref in bp'] += 1
                elif r_letter == 'N':
                    ref_uniq_array[r] = 'X'
                    self.stats['Ref N to query bp'] += 1
                else:  # No N's involved
                    if q_letter == 'X':
                        self.stats['Ref unique bp'] += 1
                    elif r_letter == 'X':
                        self.stats['Query unique bp'] += 1
                    else:
                        self.stats['Aligned Variance in bp'] += 1
            q += 1
            r += 1

        return query_uniq_array, ref_uniq_array

    def create_alignment_from_relevant_chains(self, ref_chr, relevant_chains=None):
        """In panTro4ToHg38.over.chain there are ZERO chains that have a negative strand on the reference 'tStrand'.
        I think it's a rule that you always flip the query strand instead."""
        # This assumes the chains have been sorted by score, so the highest score is the matching query_chr
        if relevant_chains is None:
            relevant_chains = [chain for chain in self.chain_list if match(chain.tName, ref_chr)]
        previous = None
        for chain in relevant_chains:
            is_master_alignment = previous is None
            self.mash_fasta_and_chain_together(chain, is_master_alignment)  # first chain is the master alignment
            previous = chain


    def setup_for_reference_chromosome(self, ref_chr, ending=''):
        ending = ref_chr + ending + \
            '__squished' * self.squish_gaps + \
            '__separate_translocations' * self.separate_translocations + \
            '__translocations' * self.show_translocations_only + \
            '__aligned_only' * self.aligned_only
        self.output_folder = make_output_dir_with_suffix(self.output_prefix, ending)
        names = {'ref': ref_chr + '_%s.fa' % first_word(self.ref_source),
                 'query': '%s_to_%s_%s.fa' % (first_word(self.query_source), first_word(self.ref_source), ref_chr)
                 }  # for collecting all the files names in a modifiable way
        # Reset values from previous iteration
        self.ref_sequence = pluck_contig(ref_chr, self.ref_source)  # only need the reference chromosome read, skip the others
        self.query_sequence = ''
        self.query_seq_gapped = array('u', '')
        self.ref_seq_gapped = array('u', '')
        self.output_fastas = []
        self.alignment = blist.blist()  # Alignment is specific to the chromosome
        self.stats = DefaultOrderedDict(lambda: 0)

        return names, ref_chr


    def _parse_chromosome_in_chain(self, chromosome_name) -> Batch:
        print("=== Begin ChainParser Unique Alignment ===")
        names, ref_chr = self.setup_for_reference_chromosome(chromosome_name)
        self.create_alignment_from_relevant_chains(ref_chr)
        self.create_fasta_from_composite_alignment()

        names['ref_gapped'], names['query_gapped'] = self.write_gapped_fasta(names['ref'], names['query'])
        names['ref_unique'], names['query_unique'] = self.print_only_unique(names['query_gapped'], names['ref_gapped'])
        # sys.exit(0)
        # NOTE: Order of these appends DOES matter!
        self.output_fastas.append(names['ref_gapped'])
        self.output_fastas.append(names['ref_unique'])
        self.output_fastas.append(names['query_unique'])
        self.output_fastas.append(names['query_gapped'])
        print("Finished creating gapped fasta files", names['ref'], names['query'])

        if True:  #self.trial_run:  # these files are never used in the viz
            del names['ref']
            del names['query']
        batch = Batch(chromosome_name, self.output_fastas, self.output_folder)
        self.write_stats_file()
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

