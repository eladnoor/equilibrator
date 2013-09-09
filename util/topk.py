"""A data structure that keeps only the top K elements it sees."""

class TopK(object):
    """Keeps the K top items."""
    
    def __init__(self, max):
        """Construction.
        
        Args:
            max: the maximum number of items to keep (k).
        """
        self.items = []
        self.max = max
    
    def _SmallestIndex(self):
        """Returns the index of the smallest element."""
        smallest = None
        smallest_index = None
        for i, item in enumerate(self.items):
            if not smallest or item < smallest:
                smallest = item
                smallest_index = i
                
        return smallest_index
        
    def GetSorted(self, key=None):
        """Return top items as a sorted list using the key function."""
        return sorted(self.items, key=key, reverse=True)

    def MaybeAdd(self, elt):
        """Potentially add elt if it's big enough.
        
        Args:
            elt: the element to add. must implement comparison.
        """
        if len(self.items) < self.max:
            self.items.append(elt)
            return
        
        # TODO(elad): cache the smallest index.
        smallest_index = self._SmallestIndex()
        if self.items[smallest_index] < elt:
            self.items[smallest_index] = elt
