"""This module deals with classes that describe a range from begin to end.
Span can have sections in the middle removed, creating two or less new Spans.
This is used by UniqueOnlyChainParser to track which parts of the file are untouched.
AlignedSpans use a pair of Span objects to track the coordinate frames of the
original and gapped sequence as gaps are added."""
from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes
from FluentDNA import gap_char

class Span(object):
    """ Span can have sections in the middle removed, creating two or less new Spans.
    This is used by UniqueOnlyChainParser to track which parts of the file are untouched."""
    def __init__(self, begin, end, contig_name=None, strand='+', zero_ok=True):
        assert zero_ok or begin != end, "%s %s are the same means zero length" % (begin, end)
        if not (begin >= 0 and end >= 0):
            raise ValueError("No negative indices! %i to %i" % (begin, end))
        self.begin = begin
        self.end = end
        self.contig_name = contig_name
        self.strand = strand
        assert self.strand in '+-'


    def __lt__(self, other_int):
        return self.begin < other_int


    def __contains__(self, index):
        if isinstance(index, Span):
            return self.overlaps(index)
        return self.begin <= index < self.end

    #def __eq__(self, other):
    #    return self.begin == other.begin and self.end == other.end and self.contig_name == other.contig_name and self.strand == other.strand

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
            if not remove_this.size():
                return None, self
            if self.size():
                raise IndexError("Remove_this doesn't overlap self at all %s %s" % (str(remove_this), str(self)))
            else:  # self has no size, just throw it away
                return None, None

        first = Span(self.begin, remove_this.begin, self.contig_name, self.strand)
        second = Span(remove_this.end, self.end, self.contig_name, self.strand)

        if remove_this.begin <= self.begin and remove_this.end >= self.end:  # delete the whole thing
            return None, None
        if remove_this.begin < self.begin < remove_this.end:  # overlaps start
            return None, second
        if remove_this.end >= self.end > remove_this.begin:  # overlaps ending
            return first, None

        return first, second  # happy path


    def sample(self, sequence, error_gap_okay=True):
        try:
            return sequence[self.begin: self.end]
        except IndexError as e:
            if error_gap_okay:
                return gap_char * len(self)
            else:
                raise e

    def set_of_points(self):
        return set(range(self.begin, self.end))



class AlignedSpans(object):
    def __init__(self, ref_span, query_span, ref_tail_size, query_tail_size,
                 is_master_chain=True, is_first_entry=False, is_hidden=False):
        """ref_span or query_span can be None to indicate an unaligned area."""
        assert ref_span.end - ref_span.begin == query_span.end - query_span.begin, "The size of the spans should be the same"
        assert query_tail_size >= 0 and ref_tail_size >= 0, "Bad tail sizes %i and %i" % (ref_tail_size, query_tail_size)
        if not (query_span.begin >= 0 and query_span.end >= 0):
            raise ValueError("No negative indices! %i to %i" % (query_span.begin, query_span.end))
        self.ref = ref_span
        self.query = query_span
        self.ref_tail_size = ref_tail_size
        self.query_tail_size = query_tail_size
        self.is_master_chain = is_master_chain
        self.is_first_entry = is_first_entry
        self.is_hidden = is_hidden


    def ref_unique_span(self):
        return Span(self.ref.end, self.ref.end + self.ref_tail_size, self.ref.contig_name, self.ref.strand)

    def query_unique_span(self):
        return Span(self.query.end, self.query.end + self.query_tail_size, self.query.contig_name, self.query.strand)


    def __lt__(self, other):
        """Useful for putting Spans in a sorted list"""
        if isinstance(other, AlignedSpans):  # TODO: more cases for None
            return self.ref.begin < other.ref.begin
        elif isinstance(other, Span):
            return self.ref.begin < other.begin
        else:
            return self.ref.begin < other

    def query_less_than(self, other):
        """Same as __lt__ but for query sequence search"""
        return self.query.begin < other.query.begin


    def __repr__(self):
        return "(%s) -> (%s)" % (str(self.ref), str(self.query))


    def align_ref_unique(self, new_alignment):
        assert isinstance(new_alignment, AlignedSpans), "This method is meant for AlignedPairs, not Spans"
        my_tail, your_tail = self.ref_unique_span().remove_from_range(new_alignment.ref)
        size = my_tail.size() if my_tail is not None else 0
        new_me = AlignedSpans(self.ref, self.query, ref_tail_size=size, query_tail_size=self.query_tail_size,
                              is_master_chain=self.is_master_chain, is_first_entry=self.is_first_entry)
        your_tail_size = your_tail.size() if your_tail is not None else 0
        try:
            you = AlignedSpans(new_alignment.ref, new_alignment.query,
                               ref_tail_size=new_alignment.ref_tail_size + your_tail_size,
                               query_tail_size=new_alignment.query_tail_size, is_master_chain=False,
                               is_first_entry=new_alignment.is_first_entry)
            return new_me, you
        except AssertionError as e:
            print(e)
            return new_me, None


    def align_query_unique(self, new_alignment):
        assert isinstance(new_alignment, AlignedSpans), "This method is meant for AlignedPairs, not Spans"
        my_tail, your_tail = self.query_unique_span().remove_from_range(new_alignment.query)
        new_me = AlignedSpans(self.ref, self.query, ref_tail_size=self.ref_tail_size,
                              query_tail_size=my_tail.size(), is_master_chain=self.is_master_chain,
                              is_first_entry=self.is_first_entry)
        you = AlignedSpans(new_alignment.ref, new_alignment.query, ref_tail_size=new_alignment.ref_tail_size,
                           query_tail_size=new_alignment.query_tail_size + your_tail.size(),
                           is_master_chain=False, is_first_entry=new_alignment.is_first_entry)
        return new_me, you


    def remove_old_query_copy(self, new_alignment):# -> tuple:
        """Each AlignedSpan also contains the record of the unaligned region following it.  In the case where
        a match has been found elsewhere in the reference, the visual representation of the sequence is moved
        to that new location based on the reference location.  That leaves behind an old, obsolete record of
        and "unalignable" region that needs to be deleted.  This function removes the old unaligned record
        without updating the query start position of the next AlignedSpan.  The end result is that the query
        sequence will be skipped over in its "native" position and must be represented in its new aligned
        location elsewhere."""
        assert isinstance(new_alignment, AlignedSpans), "This method is meant for AlignedPairs, not Spans"
        if new_alignment.query.begin in self.query_unique_span():
            my_tail, your_tail = self.query_unique_span().remove_from_range(new_alignment.query)
            new_me = AlignedSpans(self.ref, self.query, ref_tail_size=0, query_tail_size=my_tail.size(),
                                  is_master_chain=self.is_master_chain, is_first_entry=self.is_first_entry)
            progress = my_tail.size()
            you = AlignedSpans(Span(self.ref.end, self.ref.end, self.ref.contig_name, self.ref.strand),
                               Span(self.query.end + progress, self.query.end + progress,
                                    self.query.contig_name, self.query.strand), ref_tail_size=0,
                               query_tail_size=new_alignment.query.size())
            you.is_hidden = True  # This will appear as white space in the alignment
            progress += new_alignment.query.size()
            visible_tail = AlignedSpans(
                Span(self.ref.end, self.ref.end, self.ref.contig_name, self.ref.strand),
                Span(self.query.end + progress, self.query.end + progress, self.query.contig_name,
                     self.query.strand), ref_tail_size=self.ref_tail_size, query_tail_size=your_tail.size())
            return new_me, you, visible_tail
        else:
            raise IndexError(str(new_alignment.query) + " not in " + str(self.query_unique_span()))

    # def split(self, original_index):
    #     o1, o2 = self.ref_span.split(original_index)
    #     difference = original_index - self.ref_span.begin
    #     g1, g2 = self.query_span.begin + difference  # convert to gapped coordinates first
    #     first = AlignedSpan(self.ref_contig, o1.begin, o1.end)
    #     second = AlignedSpan(self.ref_contig, o2.begin, o2.end)
    #     return first, second


def alignment_chopping_index(all_alignments, new_alignment):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    lo = 0
    hi = len(all_alignments)

    while lo < hi:
        mid = (lo + hi) // 2
        if all_alignments[mid] < new_alignment:
            lo = mid + 1
        else:
            hi = mid
    return lo