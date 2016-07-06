import os


def read_seq_to_memory(file_a, file_b):
    # read contents into string for ease
    with open(file_a, 'r') as v38:
        with open(file_b, 'r') as hg19:
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


def mash_fasta_and_chain_together(query_seq, ref_seq, filename_a, filename_b):
    filename = 'hg38ToHg19.over.chain' + '__sample.chain'
    reference_is_backwards = False
    query_pointer = 0
    # ref_pointer = 0
    # ref_line_remainder = 70
    query_line_remainder = 70

    print(query_seq[:100])
    print(ref_seq[:100])
    with open(os.path.splitext(filename_a)[0] + '_gapped.fa', 'w') as query_file:
        with open(os.path.splitext(filename_b)[0] + '_gapped.fa', 'w') as ref_file:
            print('Opening', filename)
            with open(filename, 'r') as chainfile:
                lines = chainfile.readlines()
                for line_number, line in enumerate(lines):
                    if line.startswith('chain'):
                        chain, score, tName, tSize, tStrand, tStart, \
                            tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
                        query_pointer = int(qStart)
                        query_line_remainder = write_fasta_lines(query_file, 'X' * int(tStart), query_line_remainder)
                    else:
                        pieces = line.split()
                        if len(pieces) == 3:
                            size, gap_reference, gap_query = [int(x) for x in pieces]
                            if reference_is_backwards:
                                gap_reference, gap_query = gap_query, gap_reference
                            # ref_snippet = ref_seq[ref_pointer : ref_pointer + size] #+ 'X' * gap_reference
                            # ref_pointer += size #+ gap_query
                            # ref_line_remainder = write_fasta_lines(ref_file, ref_snippet, ref_line_remainder)

                            query_snippet = query_seq[query_pointer : query_pointer + size] + 'X' * gap_query
                            query_pointer += size + gap_reference
                            query_line_remainder = write_fasta_lines(query_file, query_snippet, query_line_remainder)

                        elif len(pieces) == 1:
                            # ref_line_remainder = write_fasta_lines(ref_file, ref_seq[ref_pointer: ], ref_line_remainder)
                            query_line_remainder = write_fasta_lines(query_file, query_seq[query_pointer:], query_line_remainder)


if __name__ == '__main__':
    files = ['chr9_hg38.fa', 'chr9_hg19.fa']
    v38_sequence, hg19_sequence = read_seq_to_memory(*files)
    mash_fasta_and_chain_together(v38_sequence, hg19_sequence, *files)




# def mash_fasta_and_chain_together(query_seq, gapped_filename):
#     filename = 'hg19ToHg38.over.chain' + '__sample.chain'
#     reference_is_backwards = False
#     query_pointer = 0
#     # ref_pointer = 0
#     # ref_line_remainder = 70
#     query_line_remainder = 70
#     total = 0
#     lines_processed = 0
#
#     print(query_seq[:100])
#     # with open('chr9_gapped_v38.fa', 'w') as v38_file:
#     with open(gapped_filename, 'w') as query_file:
#             print('Opening', filename)
#             with open(filename, 'r') as chainfile:
#                 lines = chainfile.readlines()
#                 for line_number, line in enumerate(lines):
#                     if line.startswith('chain'):
#                         chain, score, tName, tSize, tStrand, tStart, \
#                             tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
#                         query_pointer = int(qStart)
#                         query_line_remainder = write_fasta_lines(query_file, 'X' * int(tStart), query_line_remainder)
#                     else:
#                         pieces = line.split()
#                         if len(pieces) == 3:
#                             size, gap_reference, gap_query = [int(x) for x in pieces]
#                             if reference_is_backwards:
#                                 gap_reference, gap_query = gap_query, gap_reference
#                             # print(size, gap_reference, gap_query)
#                             snippet = query_seq[query_pointer : query_pointer + size] + 'X' * gap_query
#                             query_pointer += size + gap_reference
#                             query_line_remainder = write_fasta_lines(query_file, snippet, query_line_remainder)
#                             total += len(snippet)
#                             lines_processed += 1
#                         elif len(pieces) == 1:
#                             # ref_line_remainder = write_fasta_lines(query_file, hg19_seq[ref_pointer: ], ref_line_remainder)
#                             query_line_remainder = write_fasta_lines(query_file, query_seq[query_pointer: ], query_line_remainder)
#     print("total", total)
#     print("processed", lines_processed, len(lines), lines_processed - len(lines))