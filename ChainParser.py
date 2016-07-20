import os
from collections import defaultdict


def chunks(seq, size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


class ChainParser:
    def __init__(self):
        self.width_remaining = defaultdict(lambda: 70)
        self.ref_sequence = ''
        self.query_sequence = ''
        self.ref_seq_gapped = ''
        self.query_seq_gapped = ''

    def read_seq_to_memory(self, query_file_name, ref_file_name):
        # read contents into string for ease
        with open(query_file_name, 'r') as v38:
            with open(ref_file_name, 'r') as hg19:
                v38.readline()
                hg19.readline()  # headers
                self.query_sequence = ''.join(v38.read().splitlines())
                self.ref_sequence = ''.join(hg19.read().splitlines())


    def write_fasta_lines(self, filestream, seq):
        remaining = self.width_remaining[filestream]
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


    def mash_fasta_and_chain_together(self, chain_name, filename_a, filename_b):
        reference_is_backwards = False
        ref_collection = []
        query_collection = []
        query_pointer = 0
        ref_pointer = 0
        output_length = 0
        printed = False

        print(self.query_sequence[:100])
        print(self.ref_sequence[:100])
        gapped_fasta_query = os.path.splitext(filename_a)[0] + '_gapped.fa'
        gapped_fasta_ref = os.path.splitext(filename_b)[0] + '_gapped.fa'
        with open(gapped_fasta_query, 'w') as query_file:
            with open(gapped_fasta_ref, 'w') as ref_file:
                print('Opening', chain_name)
                with open(chain_name, 'r') as chainfile:
                    lines = chainfile.readlines()
                    for line_number, line in enumerate(lines):
                        if line.startswith('chain'):
                            chain, score, tName, tSize, tStrand, tStart, \
                                tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
                            # TODO: Show starting sequence, but have it aligned
                            query_pointer = int(tStart)  # this is correct, don't switch these around
                            ref_pointer = int(qStart)
                            # query_line_remainder = write_fasta_lines(query_file, 'X' * int(tStart), query_line_remainder)
                        else:
                            pieces = line.split()
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

                            elif len(pieces) == 1:
                                ref_collection.append(self.ref_sequence[ref_pointer:])
                                query_collection.append(self.query_sequence[query_pointer:])
                self.ref_seq_gapped = ''.join(ref_collection)
                self.query_seq_gapped = ''.join(query_collection)
                self.write_fasta_lines(ref_file, self.ref_seq_gapped)
                self.write_fasta_lines(query_file, self.query_seq_gapped)

        return gapped_fasta_ref, gapped_fasta_query


    def print_only_unique(self, ref_gapped_name, query_gapped_name):
        from array import array

        ref_unique = ref_gapped_name[:-10] + '_unique.fa'
        query_unique = query_gapped_name[:-10] + 'chimp20_unique.fa'
        overlap_size = min(len(self.ref_seq_gapped), len(self.query_seq_gapped))
        ref_array = array('u', self.ref_seq_gapped)
        que_array = array('u', self.query_seq_gapped)
        print("Done allocating unique array")
        for i in range(overlap_size):
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


    def main(self):
        import DDV
        # query_name, ref_name = 'chr9_hg19.fa', 'chr9_hg38.fa'
        # chain = 'hg19ToHg38.over.chain' + '__sample.chain'
        query_name, ref_name = 'chr20_panTro4.fa', 'chr20_hg38.fa'
        chain = 'panTro4ToHg38.over.chain__chr20_sample.chain'
        self.read_seq_to_memory(query_name, ref_name)
        ref_gapped_name, query_gapped_name = self.mash_fasta_and_chain_together(chain, query_name, ref_name)
        ref_unique_name, query_unique_name = self.print_only_unique(ref_gapped_name, query_gapped_name)

        print("Finished creating gapped fasta files", query_name, ref_name)

        DDV.DDV_main(['DDV',
                      query_gapped_name,
                      '.\\bin\\Release\\output\\dnadata\\',
                      'Parallel_Chr20_PanTro4_and_Hg38_unique_No_bg',
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
                    assert len(pieces) == 3 or len(pieces) == 1
                    chain.append(line)
    return chain


if __name__ == '__main__':
    ChainParser().main()
