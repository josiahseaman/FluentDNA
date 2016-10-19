"""This module deals with classes that describe a range from begin to end.
Span can have sections in the middle removed, creating two or less new Spans.
This is used by UniqueOnlyChainParser to track which parts of the file are untouched.
AlignedSpans use a pair of Span objects to track the coordinate frames of the
original and gapped sequence as gaps are added."""



class Span:
    """ Span can have sections in the middle removed, creating two or less new Spans.
    This is used by UniqueOnlyChainParser to track which parts of the file are untouched."""
    def __init__(self, begin, end, contig_name=None, strand='+'):
        self.begin = begin
        self.end = end
        self.contig_name = contig_name
        self.strand = strand


    def __lt__(self, other_int):
        return self.begin < other_int


    def __contains__(self, index):
        return self.begin <= index < self.end


    def __repr__(self):
        return ">%s:%i-%i" % (self.contig_name, self.begin, self.end)


    def __len__(self):
        return self.size()


    def size(self):
        return self.end - self.begin


    def overlaps(self, other):
        boundaries_check = other.begin in self or other.end - 1 in self
        is_superset = self.begin in other
        # shared_start = other.begin == self.begin
        # right_before = other.end == self.begin and other.begin != other.end
        # begin_or_end_on_wrong_side = other.end < self.begin or other.begin >= self.end
        # a = shared_start or not (begin_or_end_on_wrong_side or right_before)
        return boundaries_check or is_superset


    def split(self, split_index):
        """Splits the Span so that split_index is the first index of the second Span.
        The second span starts at split_index.  The first valid split point is begin + 1"""
        assert isinstance(self, Span), "First argument should be a Span"
        if split_index in self and split_index != self.begin + 1:
                return Span(self.begin, split_index, self.contig_name, self.strand), \
                       Span(split_index, self.end, self.contig_name, self.strand)
        raise ValueError("split_index %i is not in Span %s" % (split_index, str(self)))



class AlignedPair:
    def __init__(self, ref_span, query_span):
        """ref_span or query_span can be None to indicate an unaligned area."""
        if ref_span is not None and query_span is not None:
            assert ref_span.end - ref_span.begin == query_span.end - query_span.begin, "The size of the spans should be the same"
        self.ref = ref_span
        self.query = query_span


    def __lt__(self, other):
        """Useful for putting Spans in a sorted list"""
        if isinstance(other, AlignedPair):  # TODO: more cases for None
            return self.ref.begin < other.ref.begin
        else:
            return self.ref.begin < other

    # def split(self, original_index):
    #     o1, o2 = self.ref_span.split(original_index)
    #     difference = original_index - self.ref_span.begin
    #     g1, g2 = self.query_span.begin + difference  # convert to gapped coordinates first
    #     first = AlignedSpan(self.ref_contig, o1.begin, o1.end)
    #     second = AlignedSpan(self.ref_contig, o2.begin, o2.end)
    #     return first, second



def remove_from_range(original, remove_this):
    """Original is a range defined by (start, end).  Remove a middle range 'remove_this'
    with a (start, end) and you end up with a pair of two ranges on either side of the removal.
    Special casing for the removal overlapping the beginning or end."""
    assert isinstance(original, Span) and isinstance(remove_this, Span)

    # doesn't even overlap
    if not original.overlaps(remove_this):
        raise IndexError("Remove_this doesn't overlap original at all %s %s" % (str(remove_this), str(original)))

    first = Span(original.begin, remove_this.begin, original.contig_name, original.strand)
    second = Span(remove_this.end, original.end, original.contig_name, original.strand)

    if remove_this.begin <= original.begin and remove_this.end >= original.end:  # delete the whole thing
        return None, None
    if remove_this.begin <= original.begin < remove_this.end:  # overlaps start
        return None, second
    if remove_this.end >= original.end > remove_this.begin:  # overlaps ending
        return first, None

    return first, second  # happy path