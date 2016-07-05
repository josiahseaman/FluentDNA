from textwrap import wrap


def read_seq_to_memory():
    # read contents into string for ease
    with open('ref_GRCH38_chr9.fa', 'r') as v38:
        with open('ref_GRCH37_chr9.fa', 'r') as hg19:
            v38.readline()
            hg19.readline()  # headers
            v38_seq = ''.join(v38.read().splitlines())
            hg19_seq = ''.join(hg19.read().splitlines())
            return v38_seq, hg19_seq


def write_fasta_lines(filestream, seq, width_remaining):
    seq, remainder = seq[width_remaining:], seq[:width_remaining]
    filestream.write(remainder + '\n')
    lines = wrap(seq, width=70)
    filestream.writelines([x + '\n' for x in lines[:-1]])  # don't include last line
    left_over_size = 70 - len(lines[-1])
    filestream.write(lines[-1])  # no newline
    return left_over_size


def mash_fasta_and_chain_together(v38_seq, hg19_seq):
    filename = 'hg38ToHg19.over.chain' + '__sample.chain'
    ref_pointer, query_pointer, ref_line_remainder, query_line_remainder = 0, 0, 70, 70
    print(v38_seq[:100])
    print(hg19_seq[:100])
    with open('chr9_gapped_v38.fa', 'w') as v38_file:
        with open('chr9_gapped_hg19.fa', 'w') as hg19_file:
            print('Opening', filename)
            with open(filename, 'r') as chainfile:
                for line in chainfile.readlines():
                    if line.startswith('chain'):
                        chain, score, tName, tSize, tStrand, tStart, \
                        tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
                        ref_pointer = int(tStart)
                        query_pointer = int(qStart)
                    else:
                        pieces = line.split()
                        if len(pieces) == 3:
                            size, gap_reference, gap_query = [int(x) for x in pieces]
                            ref = v38_seq[ref_pointer: ref_pointer + size] + 'X' * gap_reference
                            ref_pointer += size
                            v38_file.write( ref + '\n')

                            query = hg19_seq[query_pointer: query_pointer + size] + 'X' * gap_query
                            query_pointer += size
                            hg19_file.write( query + '\n')

                        elif len(pieces) == 1:
                            v38_file.write( v38_seq[ref_pointer: ref_pointer + int(pieces[0])])
                            hg19_file.write( hg19_seq[query_pointer: query_pointer + int(pieces[0])])

if __name__ == '__main__':
    v38_seq, hg19_seq = read_seq_to_memory()
    mash_fasta_and_chain_together(v38_seq, hg19_seq)
