import os
from collections import defaultdict


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
    printing = False
    with open(genome_source, 'r') as genome:
        for line in genome.readlines():
            line = line.rstrip()
            if printing:
                seq_collection.append(line.upper())
            if line.startswith('>'):
                if line == chromosome_name:
                    printing = True
                    print("Found", line)
                    contig_header = line
                elif printing:
                    break  # we've collected all sequence and reached the beginning of the next contig

    return ''.join(seq_collection), contig_header


class ChainParser:
    def __init__(self):
        self.width_remaining = defaultdict(lambda: 70)
        self.ref_sequence = ''
        self.query_sequence = ''
        self.ref_seq_gapped = ''
        self.query_seq_gapped = ''


    def read_seq_to_memory(self, chromosome_name, query_source, ref_source, query_file_name, ref_file_name):
        self.query_sequence, contig_header = find_contig(chromosome_name, query_source)
        self.write_fasta_lines(query_file_name, self.query_sequence)

        self.ref_sequence, contig_header = find_contig(chromosome_name, ref_source)
        self.write_fasta_lines(ref_file_name, self.ref_sequence)


    def write_fasta_lines(self, filestream, seq):
        if isinstance(filestream, str):  # I'm actually given a file name and have to open it myself
            file_name = filestream
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

    def mash_fasta_and_chain_together(self, chain, filename_a, filename_b):
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
                # query_line_remainder = write_fasta_lines(query_file, 'X' * int(tStart), query_line_remainder)
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
        gapped_fasta_query = os.path.splitext(filename_a)[0] + '_gapped.fa'
        gapped_fasta_ref = os.path.splitext(filename_b)[0] + '_gapped.fa'
        with open(gapped_fasta_query, 'w') as query_file:
            with open(gapped_fasta_ref, 'w') as ref_file:
                self.ref_seq_gapped = ''.join(ref_collection)
                self.query_seq_gapped = ''.join(query_collection)
                self.write_fasta_lines(ref_file, self.ref_seq_gapped)
                self.write_fasta_lines(query_file, self.query_seq_gapped)
        return gapped_fasta_query, gapped_fasta_ref


    def print_only_unique(self, ref_gapped_name, query_gapped_name):
        from array import array

        ref_unique = ref_gapped_name[:-10] + '_unique.fa'
        query_unique = query_gapped_name[:-10] + '_unique.fa'
        ref_array = array('u', self.ref_seq_gapped)
        que_array = array('u', self.query_seq_gapped)
        print("Done allocating unique array")
        for i in range(min(len(self.ref_seq_gapped), len(self.query_seq_gapped))):
            # only overlapping section
            if self.ref_seq_gapped[i] == self.query_seq_gapped[i]:
                ref_array[i] = 'X'
                que_array[i] = 'X'
            else:  # Not equal
                if self.ref_seq_gapped[i] == 'N':
                    ref_array[i] = 'X'
                if self.query_seq_gapped[i] == 'N':
                    que_array[i] = 'X'

        # TODO: handle dangling on either side
        with open(ref_unique, 'w') as ref_filestream:
            with open(query_unique, 'w') as query_filestream:
                self.write_fasta_lines(ref_filestream, ''.join(ref_array))
                self.write_fasta_lines(query_filestream, ''.join(que_array))

        return ref_unique, query_unique


    def main(self, chromosome_name):
        import DDV

        query_source = 'panTro4.fa'
        ref_source = 'Human v38 All Variants.fa'

        query_name, ref_name = chromosome_name + '_panTro4.fa', chromosome_name + '_hg38.fa'
        chain = sample_chain_file(chromosome_name)
        self.read_seq_to_memory(chromosome_name, query_source, ref_source, query_name, ref_name)
        ref_gapped_name, query_gapped_name = self.mash_fasta_and_chain_together(chain, query_name, ref_name)
        ref_unique_name, query_unique_name = self.print_only_unique(ref_gapped_name, query_gapped_name)

        print("Finished creating gapped fasta files", query_name, ref_name)

        DDV.DDV_main(['DDV',
                      query_gapped_name,
                      '.\\bin\\Release\\output\\dnadata\\',
                      'Parallel_' + chromosome_name + '_PanTro4_and_Hg38__upper_case',
                      query_unique_name, ref_unique_name,
                      ref_gapped_name])


def sample_chain_file(chromosome_name, filename='panTro4ToHg38.over.chain'):
    chain = []
    with open(filename, 'r') as infile:
        printing = False
        for line in infile.readlines():
            if line.startswith('chain'):
                label, score, tName, tSize, tStrand, tStart, \
                tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
                if printing:
                    break  # only output the first chain
                printing = chromosome_name in tName and chromosome_name in qName
                if printing:
                    print(line)
                    chain.append(line)  # write header
            else:
                pieces = line.split()
                if printing:
                    if len(pieces) == 3 or len(pieces) == 1:
                        chain.append(line)
                    elif line:
                        print("bad line", line)
    return chain


if __name__ == '__main__':
    ChainParser().main('chr20')
