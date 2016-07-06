
def read_seq_to_memory():
    # read contents into string for ease
    with open('chr9_hg38.fa', 'r') as v38:
        with open('chr9_hg19.fa', 'r') as hg19:
            v38.readline()
            hg19.readline()  # headers
            v38_seq = ''.join(v38.read().splitlines())
            hg19_seq = ''.join(hg19.read().splitlines())
            return v38_seq, hg19_seq


def chunks(seq, size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def write_fasta_lines(filestream, seq, width_remaining):
    try:
        left_over_size = 70
        remainder = seq[:width_remaining]
        seq = seq[width_remaining:]
        filestream.write(remainder + '\n')
        for line in chunks(seq, 70):
            if len(line) == 70:
                filestream.write(line + '\n')
            else:
                left_over_size = 70 - len(line)
                filestream.write(line)  # no newline
        return left_over_size
    except Exception as e:
        print(e)
        return 70


def mash_fasta_and_chain_together(v38_seq, hg19_seq):
    filename = 'hg38ToHg19.over.chain' + '__sample.chain'
    query_pointer = 0
    # ref_pointer = 0
    # ref_line_remainder = 70
    query_line_remainder = 70

    print(v38_seq[:100])
    print(hg19_seq[:100])
    with open('chr9_gapped_v38.fa', 'w') as v38_file:
        # with open('chr9_gapped_hg19.fa', 'w') as hg19_file:
            print('Opening', filename)
            with open(filename, 'r') as chainfile:
                lines = chainfile.readlines()
                for line_number, line in enumerate(lines):
                    if line.startswith('chain'):
                        chain, score, tName, tSize, tStrand, tStart, \
                        tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
                        query_pointer = int(qStart)
                        query_line_remainder = write_fasta_lines(v38_file, 'X' * int(tStart), query_line_remainder)
                    else:
                        pieces = line.split()
                        if len(pieces) == 3:
                            size, gap_reference, gap_query = [int(x) for x in pieces]
                            # ref_seq = hg19_seq[ref_pointer : ref_pointer + size] #+ 'X' * gap_reference
                            # ref_pointer += size #+ gap_query
                            # ref_line_remainder = write_fasta_lines(hg19_file, ref_seq, ref_line_remainder)
                            query_seq = v38_seq[query_pointer : query_pointer + size] + 'X' * gap_query
                            query_pointer += size + gap_reference
                            query_line_remainder = write_fasta_lines(v38_file, query_seq, query_line_remainder)

                        elif len(pieces) == 1:
                            # ref_line_remainder = write_fasta_lines(hg19_file, hg19_seq[ref_pointer: ], ref_line_remainder)
                            query_line_remainder = write_fasta_lines(v38_file, v38_seq[query_pointer: ], query_line_remainder)


if __name__ == '__main__':
    v38_seq, hg19_seq = read_seq_to_memory()
    mash_fasta_and_chain_together(v38_seq, hg19_seq)
