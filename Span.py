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
        return ">%s:%s-%s" % (self.contig_name, '{:,}'.format(self.begin), '{:,}'.format(self.end))


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


    def remove_from_range(self, remove_this):
        """self is a range defined by (start, end).  Remove a middle range 'remove_this'
        with a (start, end) and you end up with a pair of two ranges on either side of the removal.
        Special casing for the removal overlapping the beginning or end."""
        assert isinstance(self, Span) and isinstance(remove_this, Span)

        # doesn't even overlap
        if not self.overlaps(remove_this):
            raise IndexError("Remove_this doesn't overlap self at all %s %s" % (str(remove_this), str(self)))

        first = Span(self.begin, remove_this.begin, self.contig_name, self.strand)
        second = Span(remove_this.end, self.end, self.contig_name, self.strand)

        if remove_this.begin <= self.begin and remove_this.end >= self.end:  # delete the whole thing
            return None, None
        if remove_this.begin <= self.begin < remove_this.end:  # overlaps start
            return None, second
        if remove_this.end >= self.end > remove_this.begin:  # overlaps ending
            return first, None

        return first, second  # happy path



class AlignedSpans:
    def __init__(self, ref_span, query_span, query_tail_size, ref_tail_size):
        """ref_span or query_span can be None to indicate an unaligned area."""
        assert ref_span.end - ref_span.begin == query_span.end - query_span.begin, "The size of the spans should be the same"
        self.ref = ref_span
        self.query = query_span
        self.ref_tail_size = ref_tail_size
        self.query_tail_size = query_tail_size

    def ref_unique_span(self):
        return Span(self.ref.end, self.ref.end + self.ref_tail_size, self.ref.contig_name, self.ref.strand)

    def query_unique_span(self):
        return Span(self.query.end, self.query.end + self.query_tail_size, self.query.contig_name, self.query.strand)

    def __lt__(self, other):
        """Useful for putting Spans in a sorted list"""
        if isinstance(other, AlignedSpans):  # TODO: more cases for None
            return self.ref.begin < other.ref.begin
        else:
            return self.ref.begin < other


    def __repr__(self):
        return "(%s) -> (%s)" % (str(self.ref), str(self.query))


    def align_ref_unique(self, new_alignment):
        assert isinstance(new_alignment, AlignedSpans), "This method is meant for AlignedPairs, not Spans"
        my_tail, your_tail = self.ref_unique_span().remove_from_range(new_alignment.ref)
        size = my_tail.size() if my_tail is not None else 0
        new_me = AlignedSpans(self.ref, self.query, query_tail_size=self.query_tail_size, ref_tail_size=size)
        your_tail_size = your_tail.size() if your_tail is not None else 0
        you = AlignedSpans(new_alignment.ref, new_alignment.query,
                           query_tail_size=new_alignment.query_tail_size,
                           ref_tail_size=new_alignment.ref_tail_size + your_tail_size)
        return new_me, you


    def align_query_unique(self, new_alignment):
        assert isinstance(new_alignment, AlignedSpans), "This method is meant for AlignedPairs, not Spans"
        my_tail, your_tail = self.query_unique_span().remove_from_range(new_alignment.query)
        new_me = AlignedSpans(self.ref, self.query, query_tail_size=my_tail.size(), ref_tail_size=self.ref_tail_size)
        you = AlignedSpans(new_alignment.ref, new_alignment.query,
                           query_tail_size=new_alignment.query_tail_size + your_tail.size(),
                           ref_tail_size=new_alignment.ref_tail_size)
        return new_me, you


    # def split(self, original_index):
    #     o1, o2 = self.ref_span.split(original_index)
    #     difference = original_index - self.ref_span.begin
    #     g1, g2 = self.query_span.begin + difference  # convert to gapped coordinates first
    #     first = AlignedSpan(self.ref_contig, o1.begin, o1.end)
    #     second = AlignedSpan(self.ref_contig, o2.begin, o2.end)
    #     return first, second


