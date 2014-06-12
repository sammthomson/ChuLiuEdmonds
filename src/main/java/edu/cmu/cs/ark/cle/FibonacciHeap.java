package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;

public class FibonacciHeap<T> {
	// the largest degree of any root list entry is guaranteed to be <= log_phi(Integer.MAX_VALUE) = 45
	private final static int MAX_DEGREE = 45;

    private Optional<Entry<T>> oMinEntry = Optional.absent();
    private int size;

	static class Entry<T> {
		private T value;
		private double priority;
		private Optional<Entry<T>> oParent = Optional.absent();
		private Optional<Entry<T>> oFirstChild = Optional.absent();
		private Entry<T> previous;
		private Entry<T> next;
		private int degree = 0;
		/** Whether this entry has had a child removed since it was
		 * added to its parent. */
		private boolean isMarked = false;

		private Entry(T value, double priority) {
			this.value = value;
			this.priority = priority;
			previous = next = this;
		}

		public static <T> Entry<T> of(T value, double priority) {
			return new Entry<T>(value, priority);
		}

		/**
		 * Cuts this from its parent and adds it to the root list, and then
		 * does the same for its parent, and so on up the tree.
		 *
		 * Runtime: O(log n)
		 */
		private void cut(Entry<T> minEntry) {
			if (!oParent.isPresent()) return;
			final Entry<T> parent = oParent.get();
			parent.degree--;
			isMarked = false;
			// update parent's `oFirstChild` pointer
			assert parent.oFirstChild.isPresent();
			if (parent.oFirstChild.get().equals(this)) {
				if (parent.degree == 0) {
					parent.oFirstChild = Optional.absent();
				} else {
					parent.oFirstChild = Optional.of(next);
				}
			}
			oParent = Optional.absent();
			unlinkFromNeighbors();
			// add to root list
			mergeLists(Optional.of(this), Optional.of(minEntry));
			if (parent.oParent.isPresent()) {
				if (parent.isMarked) {
					parent.cut(minEntry);
				} else {
					parent.isMarked = true;
				}
			}
		}

		/**
		 * Make this entry a child of the given parent entry. All linkages
		 * are updated, the degree of the parent is incremented, and
		 * `isMarked` is set to false.
		 */
		private void setParent(Entry<T> parent) {
			// remove this from its circular list
			unlinkFromNeighbors();
			oParent = Optional.of(parent);
			parent.oFirstChild = mergeLists(Optional.of(this), parent.oFirstChild);
			parent.degree++;
			isMarked = false;
		}

		private void unlinkFromNeighbors() {
			previous.next = next;
			next.previous = previous;
			previous = this;
			next = this;
		}
	}

	public static <U> FibonacciHeap<U> create() {
		return new FibonacciHeap<U>();
	}

	/**
     * Removes all elements from the heap.
     *
     * Runtime: O(1)
     */
    public void clear() {
        oMinEntry = Optional.absent();
        size = 0;
    }

    /**
     * Consolidates the trees in the heap by joining trees of equal
     * degree until there are no more trees of equal degree in the
     * root list.
     *
     * Runtime: O(log n) (amortized)
	 */
    private static <T> Optional<Entry<T>> consolidate(Entry<T> someRoot) {
		final Entry<T>[] rootsByDegree = new Entry[MAX_DEGREE];
		// For each root list entry look for others of the same degree.
		// move the larger priority root beneath the smaller priority root
		Entry<T> start = someRoot;
		Entry<T> currRoot = start;
		do {
			Entry<T> next = currRoot.next;
			int currDegree = currRoot.degree;
			Entry<T> leastRootOfCurrDegree = currRoot;
			while (rootsByDegree[currDegree] != null) {
				// Make one of the entries a child of the other.
				Entry<T> oldRootOfCurrDegree = rootsByDegree[currDegree];
				// ensure leastRootOfCurrDegree.priority <= oldRootOfCurrDegree.priority
				if (oldRootOfCurrDegree.priority < leastRootOfCurrDegree.priority) {
					final Entry<T> swap = oldRootOfCurrDegree;
					oldRootOfCurrDegree = leastRootOfCurrDegree;
					leastRootOfCurrDegree = swap;
				}
				// we're about to make `oldRootOfCurrDegree` not a root anymore, so
				// make sure we're still looking at roots
				if (start.equals(oldRootOfCurrDegree)) start = start.next;
				if (next.equals(oldRootOfCurrDegree)) next = next.next;
				// move the larger priority root beneath the smaller priority root
				oldRootOfCurrDegree.setParent(leastRootOfCurrDegree);
				// `leastRootOfCurrDegree` now has degree `currDegree` + 1
				rootsByDegree[currDegree] = null;
				currDegree++;
			}
			rootsByDegree[currDegree] = leastRootOfCurrDegree;
			currRoot = next;
		} while (currRoot != start);

		// Find the minimum priority root
		Entry<T> minRoot = start;
		for (Entry<T> root : rootsByDegree) {
			if (root != null && root.priority < minRoot.priority) {
				minRoot = root;
			}
		}
		return Optional.of(minRoot);
	}

    /**
     * Decreases the priority for an entry.
	 *
     * Runtime: O(1) (amortized)
     */
    public void decreasePriority(Entry<T> x, double priority) {
		Preconditions.checkArgument(priority <= x.priority, "Cannot increase priority");
		assert oMinEntry.isPresent();
		final Entry<T> minEntry = oMinEntry.get();
		x.priority = priority;
		if (x.oParent.isPresent() && (priority < x.oParent.get().priority)) {
			x.cut(minEntry);
        }
		if (priority < minEntry.priority) {
			oMinEntry = Optional.of(x);
		}
	}

	/**
     * Deletes `entry` from the heap.
     * The heap will be consolidated, if necessary.
     *
     * Runtime: O(log n) amortized
     */
    public void delete(Entry<T> entry) {
		assert oMinEntry.isPresent();
		final Entry<T> minEntry = oMinEntry.get();
		entry.priority = Double.POSITIVE_INFINITY;
		if (entry.oParent.isPresent()) {
			entry.cut(minEntry);
        }
		oMinEntry = Optional.of(entry);
		removeMin();
    }

    /**
     * Returns true if the heap is empty, false otherwise.
     *
     * Runtime: O(1)
     */
    public boolean isEmpty() {
        return !oMinEntry.isPresent();
    }

    /**
     * Inserts a new value element into the heap.
	 * No heap consolidation is performed.
     *
     * Runtime: O(1)
     */
    public Entry<T> insert(T value, double priority) {
		Preconditions.checkNotNull(value);
        final Entry<T> result = Entry.of(value, priority);
		oMinEntry = mergeLists(Optional.of(result), oMinEntry);
        size++;
        return result;
    }

    /**
     * Returns the entry with the minimum priority, or absent if empty.
     *
     * Runtime: O(1)
     */
    public Optional<Entry<T>> min() {
        return oMinEntry;
    }

    /**
     * Removes the smallest element from the heap. This will cause
     * the trees in the heap to be consolidated, if necessary.
     *
     * Runtime: O(log n) amortized
     *
     * @return  value object with the smallest priority.
     */
    public Optional<T> removeMin() {
		if (!oMinEntry.isPresent()) {
			return Optional.absent();
		}
		final Entry<T> minEntry = oMinEntry.get();
		// move minEntry's children to the root list
		if (minEntry.oFirstChild.isPresent()) {
			final Entry<T> minFirstChild = minEntry.oFirstChild.get();
			minFirstChild.oParent = Optional.absent();
            for (Entry<T> childOfMin = minFirstChild.next; childOfMin != minFirstChild; childOfMin = childOfMin.next) {
                childOfMin.oParent = Optional.absent();
            }
			mergeLists(oMinEntry, minEntry.oFirstChild);
        }
        // remove minEntry from root list
		final Entry<T> next = minEntry.next;
		if (minEntry.equals(next)) {
            oMinEntry = Optional.absent();
        } else {
			minEntry.unlinkFromNeighbors();
			oMinEntry = consolidate(next);
        }
        size--;
        return Optional.of(minEntry.value);
    }

    /**
     * Returns the number of elements in the heap.
     *
     * Runtime: O(1)
     */
    public int size() {
        return size;
    }

    /**
     * Joins two Fibonacci heaps into a new one. No heap consolidation is
     * performed; the two root lists are just spliced together.
     *
     * Runtime: O(1)
     */
    public static <U> FibonacciHeap<U> merge(FibonacciHeap<U> a, FibonacciHeap<U> b) {
        final FibonacciHeap<U> result = FibonacciHeap.create();
		result.oMinEntry = mergeLists(a.oMinEntry, b.oMinEntry);
		result.size = a.size + b.size;
		return result;
    }

	/**
	 * Merge two doubly-linked circular lists, given a pointer into each.
	 * Return the smaller of the two arguments.
	 */
	private static <T> Optional<Entry<T>> mergeLists(Optional<Entry<T>> oA, Optional<Entry<T>> oB) {
		if (!oA.isPresent()) return oB;
		if (!oB.isPresent()) return oA;
		final Entry<T> a = oA.get();
		final Entry<T> b = oB.get();
		// splice the two circular lists together like a Mobius strip
		final Entry<T> aOldNext = a.next;
		a.next = b.next;
		a.next.previous = a;
		b.next = aOldNext;
		b.next.previous = b;
		return a.priority < b.priority ? oA : oB;
	}
}
