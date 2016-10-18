"""This module deals with classes that describe a range from begin to end.
Span can have sections in the middle removed, creating two or less new Spans.
This is used by UniqueOnlyChainParser to track which parts of the file are untouched.
AlignedSpans use a pair of Span objects to track the coordinate frames of the
original and gapped sequence as gaps are added."""



class Span:
    """ Span can have sections in the middle removed, creating two or less new Spans.
    This is used by UniqueOnlyChainParser to track which parts of the file are untouched."""
    def __init__(self, begin, end):
        self.begin = begin
        self.end = end

    def __lt__(self, other_int):
        return self.begin < other_int

    def __contains__(self, index):
        return self.begin <= index < self.end

    def overlaps(self, other):
        boundaries_check = other.begin in self or other.end - 1 in self
        is_superset = self.begin in other
        # shared_start = other.begin == self.begin
        # right_before = other.end == self.begin and other.begin != other.end
        # begin_or_end_on_wrong_side = other.end < self.begin or other.begin >= self.end
        # a = shared_start or not (begin_or_end_on_wrong_side or right_before)
        return boundaries_check or is_superset

    def split(self, split_index):
        assert isinstance(self, Span), "First argument should be a Span"
        if split_index in self:
            return Span(self.begin, split_index + 1), Span(split_index, self.end)
        raise ValueError("split_index %i is not in Span %s" % (split_index, str(self)))



class AlignedSpan:
    def __init__(self, original_begin, original_end, gapped_begin, gapped_end):
        assert original_end - original_begin == gapped_end - gapped_begin, "The size of the spans should be the same"
        self.original = Span(original_begin, original_end)
        self.gapped = Span(gapped_begin, gapped_end)

    def __lt__(self, other):
        """Useful for putting Spans in a sorted list"""
        if isinstance(other, AlignedSpan):
            return self.original.begin < other.original.begin
        else:
            return self.original.begin < other

    def split(self, original_index):
        o1, o2 = self.original.split(original_index)
        difference = original_index - self.original.begin
        g1, g2 = self.gapped.begin + difference  # convert to gapped coordinates first
        return AlignedSpan(o1.begin, o1.end, g1.begin, g1.end), AlignedSpan(o2.begin, o2.end, g2.begin, g2.end)



def remove_from_range(original, remove_this):
    """Original is a range defined by (start, end).  Remove a middle range 'remove_this'
    with a (start, end) and you end up with a pair of two ranges on either side of the removal.
    Special casing for the removal overlapping the beginning or end."""
    assert isinstance(original, Span) and isinstance(remove_this, Span)

    # doesn't even overlap
    if not original.overlaps(remove_this):
        raise IndexError("Remove_this doesn't overlap original at all %s %s" % (str(remove_this), str(original)))

    first = Span(original.begin, remove_this.begin)
    second = Span(remove_this.end, original.end)

    if remove_this.begin <= original.begin and remove_this.end >= original.end:  # delete the whole thing
        return None, None
    if remove_this.begin <= original.begin < remove_this.end:  # overlaps start
        return None, second
    if remove_this.end >= original.end > remove_this.begin:  # overlaps ending
        return first, None

    return first, second  # happy path
