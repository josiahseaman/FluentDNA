import os
import shutil
import multiprocessing

from collections import defaultdict
from array import array


def chunks(seq, size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def find_contig(chromosome_name, genome_source):
    """Scan through a genome fasta file looking for a matching contig name.  When it find it, find_contig collects
    the sequence and returns it as a string with no cruft."""
    chromosome_name = '>' + chromosome_name
    contig_header = ''
    seq_collection = []
    headers = []
    printing = False
    with open(genome_source, 'r') as genome:
        for line in genome.readlines():
            line = line.rstrip()
            if printing:
                seq_collection.append(line.upper())  # always upper case so equality checks work
            if line.startswith('>'):
                # headers.append(line)
                if line == chromosome_name:
                    printing = True
                    print("Found", line)
                    contig_header = line
                elif printing:
                    break  # we've collected all sequence and reached the beginning of the next contig
    assert len(seq_collection), "Contig not found." + chromosome_name  # File contained these contigs:\n" + '\n'.join(headers)
    return ''.join(seq_collection), contig_header


def sample_chain_file(chromosome_name, filename):
    chain = []
    with open(filename, 'r') as infile:
        printing = False
        for line in infile.readlines():
            if line.startswith('chain'):
                label, score, tName, tSize, tStrand, tStart, \
                tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
                if printing:
                    break  # only output the first chain
                printing = chromosome_name in tName and chromosome_name in qName and '+' in tStrand and '+' in qStrand
                if printing:
                    print(line)
                    chain.append(line)  # write header
            else:
                pieces = line.split()
                if printing:
                    if len(pieces) == 3 or len(pieces) == 1:
                        chain.append(line)
    return chain


class ChainParser:
    extra_generated_fastas = []

    def __init__(self, fasta, comparison_fasta, chain_file, output_folder):
        self.width_remaining = defaultdict(lambda: 70)
        self.ref_sequence = ''
        self.query_sequence = ''
        self.ref_seq_gapped = ''
        self.query_seq_gapped = ''
        self.fasta = fasta
        self.comparison_fasta = comparison_fasta
        self.chain_file = chain_file
        self.output_folder = output_folder

    def read_seq_to_memory(self, chromosome_name, query_source, ref_source, query_file_name, ref_file_name):
        self.query_sequence, contig_header = find_contig(chromosome_name, query_source)
        self.write_fasta_lines(query_file_name, self.query_sequence)

        self.ref_sequence, contig_header = find_contig(chromosome_name, ref_source)
        self.write_fasta_lines(ref_file_name, self.ref_sequence)

    def write_fasta_lines(self, filestream, seq):
        if isinstance(filestream, str):  # I'm actually given a file name and have to open it myself
            file_name = os.path.join(self.output_folder, filestream)
            with open(file_name, 'w') as filestream:
                self.do_write(filestream, 70, seq)
                print("Wrote", file_name, len(self.query_sequence))
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

    def mash_fasta_and_chain_together(self, chain, filename_a, filename_b, is_master_alignment=False):
        reference_is_backwards = False
        ref_collection = []
        query_collection = []
        query_pointer = 0
        ref_pointer = 0
        output_length = 0
        printed = False

        print(self.query_sequence[:100])
        print(self.ref_sequence[:100])
        for chain_line in chain:
            if chain_line.startswith('chain'):
                label, score, tName, tSize, tStrand, tStart, \
                tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = chain_line.split()
                # TODO: Show starting sequence, but have it aligned
                query_pointer = int(tStart)  # this is correct, don't switch these around
                ref_pointer = int(qStart)
                if is_master_alignment:  # include the unaligned beginning of the sequence
                    longer_gap = max(ref_pointer, query_pointer)
                    ref_collection.append(self.ref_sequence[:ref_pointer] + 'X' * (longer_gap - ref_pointer))  # one of these gaps will be 0
                    query_collection.append(self.query_sequence[:query_pointer] + 'X' * (longer_gap - query_pointer))
            else:
                pieces = chain_line.split()
                if len(pieces) == 3:
                    size, gap_reference, gap_query = [int(x) for x in pieces]
                    if reference_is_backwards:
                        gap_reference, gap_query = gap_query, gap_reference

                    if not printed and output_length > 4319440:
                        printed = True
                        print("Start at", size, gap_reference, gap_query)
                        print(filename_a, query_pointer)
                        print(filename_b, ref_pointer)
                        print(output_length)

                    space_saved = max(0, min(gap_query, gap_reference))

                    ref_snippet = self.ref_sequence[ref_pointer: ref_pointer + size + gap_query] + 'X' * (gap_reference - space_saved)
                    ref_pointer += size + gap_query  # alignable and unalignable block concatenated together
                    ref_collection.append(ref_snippet)
                    output_length += len(ref_snippet)

                    query_snippet = self.query_sequence[query_pointer: query_pointer + size] + 'X' * (gap_query - space_saved)
                    query_snippet += self.query_sequence[query_pointer + size: query_pointer + size + gap_reference]
                    assert len(query_snippet) == len(ref_snippet), "You should be outputting equal length strings til the end"
                    query_pointer += size + gap_reference  # two blocks of sequence separated by gap
                    query_collection.append(query_snippet)

                elif len(pieces) == 1:  # last one: print out all remaining sequence
                    ref_collection.append(self.ref_sequence[ref_pointer:])
                    query_collection.append(self.query_sequence[query_pointer:])

        gapped_fasta_query, gapped_fasta_ref = self.write_gapped_fasta(filename_a, filename_b, query_collection, ref_collection)

        return gapped_fasta_ref, gapped_fasta_query

    def write_gapped_fasta(self, filename_a, filename_b, query_collection, ref_collection):
        gapped_fasta_query = os.path.join(self.output_folder, os.path.splitext(filename_a)[0] + '_gapped.fa')
        gapped_fasta_ref = os.path.join(self.output_folder, os.path.splitext(filename_b)[0] + '_gapped.fa')
        with open(gapped_fasta_query, 'w') as query_file:
            with open(gapped_fasta_ref, 'w') as ref_file:
                self.ref_seq_gapped = ''.join(ref_collection)
                self.query_seq_gapped = ''.join(query_collection)
                self.write_fasta_lines(ref_file, self.ref_seq_gapped)
                self.write_fasta_lines(query_file, self.query_seq_gapped)
        self.extra_generated_fastas.append(gapped_fasta_query)
        self.extra_generated_fastas.append(gapped_fasta_ref)
        return gapped_fasta_query, gapped_fasta_ref

    def print_only_unique(self, ref_gapped_name, query_gapped_name, is_master_alignment=True):
        ref_unique = ref_gapped_name[:-10] + '_unique.fa'
        query_unique = query_gapped_name[:-10] + '_unique.fa'
        ref_array = array('u', self.ref_seq_gapped)
        que_array = array('u', self.query_seq_gapped)
        print("Done allocating unique array")
        shortest_sequence = min(len(self.ref_seq_gapped), len(self.query_seq_gapped))
        for i in range(shortest_sequence):
            # only overlapping section
            if self.ref_seq_gapped[i] == self.query_seq_gapped[i]:
                ref_array[i] = 'X'
                que_array[i] = 'X'
            else:  # Not equal
                if self.ref_seq_gapped[i] == 'N':
                    ref_array[i] = 'X'
                if self.query_seq_gapped[i] == 'N':
                    que_array[i] = 'X'

        # Just to be thorough: prints aligned section (shortest_sequence) plus any dangling end sequence
        ref_unique_filepath = os.path.join(self.output_folder, ref_unique)
        query_unique_filepath = os.path.join(self.output_folder, query_unique)
        with open(ref_unique_filepath, 'w') as ref_filestream:
            self.write_fasta_lines(ref_filestream, ''.join(ref_array))
        with open(query_unique_filepath, 'w') as query_filestream:
            self.write_fasta_lines(query_filestream, ''.join(que_array))
        self.extra_generated_fastas.append(ref_unique_filepath)
        self.extra_generated_fastas.append(query_unique_filepath)

        return ref_unique, query_unique

    def _parse_chromosome_in_chain(self, chromosome_name):
        fasta = {'query_name': chromosome_name + '_' + os.path.basename(self.fasta), 'ref_name': chromosome_name + '_' + os.path.basename(self.comparison_fasta)}  # for collecting all the files names in a modifiable way

        chain = sample_chain_file(chromosome_name, self.chain_file)
        self.read_seq_to_memory(chromosome_name, self.fasta, self.comparison_fasta, fasta['query_name'], fasta['ref_name'])
        fasta['ref_gapped_name'], fasta['query_gapped_name'] = self.mash_fasta_and_chain_together(chain, fasta['query_name'], fasta['ref_name'], True)
        fasta['ref_unique_name'], fasta['query_unique_name'] = self.print_only_unique(fasta['ref_gapped_name'], fasta['query_gapped_name'])
        print("Finished creating gapped fasta files", fasta['query_name'], fasta['ref_name'])

    def parse_chain(self, chromosomes=None):  # TODO: Remove ability to not pass in chromosomes
        assert isinstance(chromosomes, list), "'Chromosomes' must be a list! A single element list is okay."

        # TODO: handle Chr2A and Chr2B separately
        for chr in chromosomes:
            self._parse_chromosome_in_chain(chr)

