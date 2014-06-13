package edu.cmu.cs.ark.cle.ds;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

import java.util.Comparator;
import java.util.List;

/**
 * A Fibonacci heap, as described by Fredman and Tarjan.
 *
 * @param <V> the type of the values stored in the heap
 * @param <P> the type of the priorities
 */
public class FibonacciHeap<V,P> {
	public final static int MAX_CAPACITY = Integer.MAX_VALUE;
	// the largest degree of any root list entry is guaranteed to be <= log_phi(MAX_CAPACITY) = 45
	private final static int MAX_DEGREE = 45;

    private Optional<Entry> oMinEntry = Optional.absent();
    private int size = 0;
	private final Comparator<? super P> comparator;

	class Entry {
		private V value;
		private P priority;
		private Optional<Entry> oParent = Optional.absent();
		private Optional<Entry> oFirstChild = Optional.absent();
		private Entry previous;
		private Entry next;
		private int degree = 0;
		/** Whether this entry has had a child removed since it was
		 * added to its parent. */
		private boolean isMarked = false;

		private Entry(V value, P priority) {
			this.value = value;
			this.priority = priority;
			previous = next = this;
		}

	}

	private FibonacciHeap(Comparator<? super P> comparator) {
		// we'll use nulls to force a node to the top when we delete it
		this.comparator = Ordering.from(comparator).nullsFirst();
	}

	public static <T,C> FibonacciHeap<T,C> create(Comparator<? super C> comparator) {
		return new FibonacciHeap<T, C>(comparator);
	}

	public static <T,C extends Comparable> FibonacciHeap<T,C> create() {
		return FibonacciHeap.create(Ordering.<C>natural());
	}

	public Entry entry(V value, P priority) {
		return new Entry(value, priority);
	}

	/**
	 * Returns the comparator used to order the elements in this
	 * queue.
	 */
	public Comparator<? super P> comparator() {
		return comparator;
	}

	/**
     * Removes all elements from the heap.
     * Runtime: O(1)
     */
    public void clear() {
        oMinEntry = Optional.absent();
        size = 0;
    }

    /**
     * Decreases the priority for an entry.
     * Runtime: O(1) (amortized)
     */
    public void decreasePriority(Entry entry, P newPriority) {
		Preconditions.checkArgument(comparator.compare(newPriority, entry.priority) <= 0,
				"Cannot increase priority");
		entry.priority = newPriority;
		assert oMinEntry.isPresent();
		final Entry minEntry = oMinEntry.get();
		if (entry.oParent.isPresent() && (comparator.compare(newPriority, entry.oParent.get().priority) < 0)) {
			cutAndMakeRoot(entry);
        }
		if (comparator.compare(newPriority, minEntry.priority) < 0) {
			oMinEntry = Optional.of(entry);
		}
	}

	/**
     * Deletes `entry` from the heap. The heap will be consolidated, if necessary.
     * Runtime: O(log n) amortized
     */
    public void delete(Entry entry) {
		entry.priority = null;
		cutAndMakeRoot(entry);
		oMinEntry = Optional.of(entry);
		poll();
    }

    /**
     * Returns true if the heap is empty, false otherwise.
     * Runtime: O(1)
     */
    public boolean isEmpty() {
        return !oMinEntry.isPresent();
    }

    /**
     * Inserts a new value element into the heap.
	 * No heap consolidation is performed.
     * Runtime: O(1)
     */
    public Entry add(V value, P priority) {
		Preconditions.checkNotNull(value);
		Preconditions.checkNotNull(priority);
		Preconditions.checkArgument(size < MAX_CAPACITY, "Maximum capacity reached.");
        final Entry result = entry(value, priority);
		// add as a root
		oMinEntry = mergeLists(Optional.of(result), oMinEntry);
        size++;
        return result;
    }

	/**
	 * Returns the entry with the minimum priority, or absent if empty.
	 * Runtime: O(1)
	 */
	public Optional<Entry> peek() {
		return oMinEntry;
	}

	/**
     * Removes the smallest element from the heap. This will cause
     * the trees in the heap to be consolidated, if necessary.
     * Runtime: O(log n) amortized
     */
    public Optional<V> poll() {
		if (!oMinEntry.isPresent()) {
			return Optional.absent();
		}
		final Entry minEntry = oMinEntry.get();
		// move minEntry's children to the root list
		if (minEntry.oFirstChild.isPresent()) {
			final Entry minFirstChild = minEntry.oFirstChild.get();
			minFirstChild.oParent = Optional.absent();
            for (Entry childOfMin = minFirstChild.next; childOfMin != minFirstChild; childOfMin = childOfMin.next) {
                childOfMin.oParent = Optional.absent();
            }
			mergeLists(oMinEntry, minEntry.oFirstChild);
        }
        // remove minEntry from root list
		final Entry next = minEntry.next;
		if (minEntry.equals(next)) {
            oMinEntry = Optional.absent();
        } else {
			unlinkFromNeighbors(minEntry);
			oMinEntry = Optional.of(consolidate(next));
        }
        size--;
        return Optional.of(minEntry.value);
    }

	/**
     * Returns the number of elements in the heap.
     * Runtime: O(1)
     */
    public int size() {
        return size;
    }

    /**
     * Joins two Fibonacci heaps into a new one. No heap consolidation is
     * performed; the two root lists are just spliced together.
     * Runtime: O(1)
     */
    public static <U,P> FibonacciHeap<U,P> merge(FibonacciHeap<U,P> a, FibonacciHeap<U,P> b) {
		Preconditions.checkArgument(a.comparator().equals(b.comparator()),
				"Heaps that use different comparators can't be merged.");
		final FibonacciHeap<U,P> result = FibonacciHeap.create(a.comparator);
		result.oMinEntry = a.mergeLists(a.oMinEntry, b.oMinEntry);
		result.size = a.size + b.size;
		return result;
    }

	private List<Entry> getCycle(Entry start) {
		final List<Entry> results = Lists.newLinkedList();
		Entry current = start;
		do {
			results.add(current);
			current = current.next;
		} while (!current.equals(start));
		return results;
	}

	/**
	 * Merge two doubly-linked circular lists, given a pointer into each.
	 * Return the smaller of the two arguments.
	 */
	private Optional<Entry> mergeLists(Optional<Entry> oA, Optional<Entry> oB) {
		if (!oA.isPresent()) return oB;
		if (!oB.isPresent()) return oA;
		final Entry a = oA.get();
		final Entry b = oB.get();
		// splice the two circular lists together like a Mobius strip
		final Entry aOldNext = a.next;
		a.next = b.next;
		a.next.previous = a;
		b.next = aOldNext;
		b.next.previous = b;
		return comparator.compare(a.priority, b.priority) < 0 ? oA : oB;
	}

	/**
	 * Cuts this entry from its parent and adds it to the root list, and then
	 * does the same for its parent, and so on up the tree.
	 * Runtime: O(log n)
	 */
	private void cutAndMakeRoot(Entry entry) {
		if (!entry.oParent.isPresent()) return;  // already a root
		final Entry parent = entry.oParent.get();
		parent.degree--;
		entry.isMarked = false;
		// update parent's `oFirstChild` pointer
		assert parent.oFirstChild.isPresent();
		if (parent.oFirstChild.get().equals(entry)) {
			if (parent.degree == 0) {
				parent.oFirstChild = Optional.absent();
			} else {
				parent.oFirstChild = Optional.of(entry.next);
			}
		}
		entry.oParent = Optional.absent();
		unlinkFromNeighbors(entry);
		// add to root list
		mergeLists(Optional.of(entry), oMinEntry);
		if (parent.oParent.isPresent()) {
			if (parent.isMarked) {
				cutAndMakeRoot(parent);
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
	private void setParent(Entry entry, Entry parent) {
		unlinkFromNeighbors(entry);
		entry.oParent = Optional.of(parent);
		parent.oFirstChild = mergeLists(Optional.of(entry), parent.oFirstChild);
		parent.degree++;
		entry.isMarked = false;
	}

	private static void unlinkFromNeighbors(FibonacciHeap.Entry entry) {
		entry.previous.next = entry.next;
		entry.next.previous = entry.previous;
		entry.previous = entry;
		entry.next = entry;
	}

	/**
	 * Consolidates the trees in the heap by joining trees of equal
	 * degree until there are no more trees of equal degree in the
	 * root list. Returns the new minimum root.
	 * Runtime: O(log n) (amortized)
	 */
	private Entry consolidate(Entry someRoot) {
		Entry minRoot = someRoot;
		// For each root list entry look for others of the same degree.
		// move the larger priority root beneath the smaller priority root
		final Object[] rootsByDegree = new Object[MAX_DEGREE];
		for (Entry currRoot : getCycle(someRoot)) {
			int degree = currRoot.degree;
			Entry leastRootOfDegree = currRoot;
			while (rootsByDegree[degree] != null) {
				@SuppressWarnings("unchecked")
				Entry oldRootOfCurrDegree = (Entry) rootsByDegree[degree];
				// move the larger priority root beneath the smaller priority root
				if (comparator.compare(leastRootOfDegree.priority, oldRootOfCurrDegree.priority) < 0) {
					setParent(oldRootOfCurrDegree, leastRootOfDegree);
				} else {
					setParent(leastRootOfDegree, oldRootOfCurrDegree);
					leastRootOfDegree = oldRootOfCurrDegree;
				}
				// `leastRootOfDegree` now has degree `degree` + 1
				rootsByDegree[degree] = null;
				degree++;
			}
			rootsByDegree[degree] = leastRootOfDegree;
			if (comparator.compare(leastRootOfDegree.priority, minRoot.priority) <= 0) {
				minRoot = leastRootOfDegree;
			}
		}
		return minRoot;
	}
}
