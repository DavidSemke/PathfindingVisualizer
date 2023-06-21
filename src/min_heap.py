# The code for the original min_heap implementation can be found at
# https://github.com/python/cpython/blob/master/Lib/heapq.py

__all__ = ['heappush', 'heappop', 'heapify', 'heapreplace', 'heappushpop', 'heapremove']

def heappush(heap, item):
    """Push item onto heap, maintaining the heap invariant."""
    heap.append(item)
    percolations = _siftdown(heap, 0, len(heap)-1)
    return percolations

def heappop(heap):
    """Pop the smallest item off the heap, maintaining the heap invariant."""
    lastelt = heap.pop()    # raises appropriate IndexError if heap is empty
    percolations = 0
    if heap:
        returnitem = heap[0]
        heap[0] = lastelt
        percolations = _siftup(heap, 0)
        return returnitem, percolations
    return lastelt, percolations


def heapremove(heap, index):
    percolates = 0

    heap[index] = heap[-1]
    heap.pop()
    if index < len(heap):
        percolates += _siftup(heap, index)
        percolates += _siftdown(heap, 0, index)

    return percolates


def heapreplace(heap, item):
    """Pop and return the current smallest value, and add the new item.
    This is more efficient than heappop() followed by heappush(), and can be
    more appropriate when using a fixed-size heap.  Note that the value
    returned may be larger than item!  That constrains reasonable uses of
    this routine unless written as part of a conditional replacement:
        if item > heap[0]:
            item = heapreplace(heap, item)
    """
    returnitem = heap[0]    # raises appropriate IndexError if heap is empty
    heap[0] = item
    percolations = _siftup(heap, 0)
    return returnitem, percolations

def heappushpop(heap, item):
    """Fast version of a heappush followed by a heappop."""
    percolations = 0
    if heap and heap[0] < item:
        item, heap[0] = heap[0], item
        percolations += _siftup(heap, 0)
    return item, percolations

def heapify(x):
    """Transform list into a heap, in-place, in O(len(x)) time."""
    percolations = 0
    n = len(x)
    # Transform bottom-up.  The largest index there's any point to looking at
    # is the largest with a child index in-range, so must have 2*i + 1 < n,
    # or i < (n-1)/2.  If n is even = 2*j, this is (2*j-1)/2 = j-1/2 so
    # j-1 is the largest, which is n//2 - 1.  If n is odd = 2*j+1, this is
    # (2*j+1-1)/2 = j so j-1 is the largest, and that's again n//2-1.
    for i in reversed(range(n//2)):
        percolations += _siftup(x, i)
    return percolations

# 'heap' is a heap at all indices >= startpos, except possibly for pos.  pos
# is the index of a leaf with a possibly out-of-order value.  Restore the
# heap invariant.
def _siftdown(heap, startpos, pos):
    newitem = heap[pos]
    # Follow the path to the root, moving parents down until finding a place
    # newitem fits.
    percolations = 0
    while pos > startpos:
        parentpos = (pos - 1) >> 1
        parent = heap[parentpos]
        percolations += 1
        if newitem < parent:
            heap[pos] = parent
            pos = parentpos
            continue
        break
    heap[pos] = newitem

    return percolations


def _siftup(heap, pos):
    endpos = len(heap)
    startpos = pos
    newitem = heap[pos]
    percolations = 0
    # Bubble up the smaller child until hitting a leaf.
    childpos = 2*pos + 1    # leftmost child position
    while childpos < endpos:
        # Set childpos to index of smaller child.
        rightpos = childpos + 1
        if rightpos < endpos and not heap[childpos] < heap[rightpos]:
            childpos = rightpos
        # Move the smaller child up.
        heap[pos] = heap[childpos]
        percolations += 1
        pos = childpos
        childpos = 2*pos + 1
    # The leaf at pos is empty now.  Put newitem there, and bubble it up
    # to its final resting place (by sifting its parents down).
    heap[pos] = newitem
    percolations += _siftdown(heap, startpos, pos)

    return percolations
