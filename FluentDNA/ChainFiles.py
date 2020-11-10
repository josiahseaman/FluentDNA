class ChainEntry(object):
    def __init__(self, size, gap_query, gap_ref):
        self.size, self.gap_query, self.gap_ref = size, gap_query, gap_ref

    def __repr__(self):
        return "(Size: %i gap_query:%i gap_ref:%i)" % (self.size, self.gap_query, self.gap_ref)



class Chain(object):
    def __init__(self, line):
        label, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, chain_id = line.split()
        self.score = int(score)
        self.tName = tName
        self.tSize = int(tSize)
        self.tStrand = tStrand
        self.tStart = int(tStart)
        self.tEnd = int(tEnd)
        self.qName = qName
        self.qSize = int(qSize)
        self.qStrand = qStrand
        self.qStart = int(qStart)
        self.qEnd = int(qEnd)
        self.chain_id = int(chain_id)

        self.entries = []


    def __str__(self):
        return 'chain %i %s %i %s %i %i %s %i %s %i %i %i' % (self.score, self.tName, self.tSize, self.tStrand, self.tStart, self.tEnd,
                                                              self.qName, self.qSize, self.qStrand, self.qStart, self.qEnd, self.chain_id)

    def add_entry(self, line):
        pieces = line.split()
        size, gap_query, gap_ref = 0, 0, 0
        if len(pieces) == 3:
            size, gap_query, gap_ref = [int(x) for x in pieces]
        elif len(pieces) == 1:
            size = int(pieces[0])
        elif len(pieces):  # non-empty
            raise ValueError("Don't know how to parse line: " + line)
        self.entries.append(ChainEntry(size, gap_query, gap_ref))


def chain_file_to_list(chain_name, extract_contigs=None):
    """Return a list of Chain objects from a liftover .chain filename"""
    all_chains = []
    new_chain = None
    skipping_current_chain = False
    with open(chain_name, 'r') as infile:
        for line in infile.readlines():
            if line.startswith('#'):
                continue  # this is a comment, no seriously, the line in the chain file is commented out
            if line.startswith('chain'):
                if new_chain is not None:
                    all_chains.append(new_chain)
                    new_chain = None
                current_chain = Chain(line)
                if extract_contigs is None \
                        or current_chain.qName in extract_contigs \
                        or current_chain.tName in extract_contigs:
                    new_chain = current_chain
                    skipping_current_chain = False
                else:
                    skipping_current_chain = True
            elif not skipping_current_chain:
                new_chain.add_entry(line)
    if new_chain is not None:
        all_chains.append(new_chain)
    all_chains.sort(key=lambda chain: -chain.score)  # biggest score first
    return all_chains


def fetch_all_chains(ref_chr, query_chr, query_strand, chain_list, ref_strand='+'):
    """Fetches all chains that match the requirements.
    any of {query_chr, query_strand, ref_chr, ref_strand} can be None which means the match to anything.
    """
    collected = []
    for chain in chain_list:
        if match(ref_chr, chain.tName) and match(query_chr, chain.qName) and match(ref_strand, chain.tStrand) and match(query_strand, chain.qStrand):
            collected.append(chain)
    print(ref_chr, query_chr, query_strand, ":", len(collected))
    return collected


def match(target, current):
    """Returns true if current satisfies the target requirements"""
    if target is None:
        return True
    first_part = target.split()[0].upper()
    scrubbed = current.split()[0].upper()
    return target.upper() == current.upper() or first_part == scrubbed or \
           first_part == scrubbed.replace(';', '').replace('=', '')
