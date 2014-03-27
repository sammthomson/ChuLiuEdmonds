package edu.cmu.cs.ark.cle;

/**
 * Union-Find Data Structure
 *
 * @author sthomson@cs.cmu.edu
 */
public class Partition {
	private final int[] parents;
	private final int[] ranks;

	/**
	 * Constructs a new partition of singletons
	 * @param n the size of your collection
	 */
	public Partition(int n) {
		parents = new int[n];
		ranks = new int[n];
		for(int i = 0; i < n; i++) {
			parents[i] = i; // each node is its own head
			ranks[i] = 0; // every node has depth 0 to start
		}
	}

	/** Find the representative for the given item */
	public int componentOf(int i) {
		if (parents[i] == i) {
			return i;
		} else {
			parents[i] = componentOf(parents[i]);
		}
		return parents[i];
	}

	/** Merges the given components and returns the representative of the new component */
	public int merge(int a, int b) {
		final int aHead = componentOf(a);
		final int bHead = componentOf(b);
		if(aHead == bHead) return aHead;
		// add the shorter tree underneath the taller tree
		final int aRank = ranks[aHead];
		final int bRank = ranks[bHead];
		if (aRank > bRank) {
			parents[bHead] = aHead;
			return aHead;
		} else if (bRank > aRank) {
			parents[aHead] = bHead;
			return bHead;
		} else {
			// whoops, the tree got taller
			parents[bHead] = aHead;
			ranks[aHead] = aRank + 1;
			return aHead;
		}
	}

	/** Determines whether the two items are in the same component or not */
	public boolean sameComponent(int a, int b) {
		return componentOf(a) == componentOf(b);
	}
}
