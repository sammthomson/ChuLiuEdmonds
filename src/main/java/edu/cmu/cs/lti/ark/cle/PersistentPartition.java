package edu.cmu.cs.lti.ark.cle;

import com.google.common.collect.Lists;

import java.util.List;

/**
 * Implementation of a persistent Union-Find Data Structure.
 * I.e 'merge(a, b)' returns a new PersistentPartition,
 * leaving the original PersistentPartition unchanged.
 *
 * Based on "A Persistent Union-Find Data Structure" Conchon & Filliaˆtre, 2007
 *
 * NOT THREADSAFE.
 * Accesses to a PersistentArrayList can potentially mutate a chain of PersistentArrayLists
 * in a thread-unsafe way during the method call.
 *
 * @author sthomson@cs.cmu.edu
 */
public class PersistentPartition {
	private PersistentArrayList<Integer> parents;
	private final PersistentArrayList<Integer> ranks;

	/**
	 * Constructs a new partition of singletons
	 * @param n the size of your collection
	 */
	public PersistentPartition(int n) {
		final List<Integer> parents = Lists.newArrayListWithExpectedSize(n);
		final List<Integer> ranks = Lists.newArrayListWithExpectedSize(n);
		for(int i = 0; i < n; i++) {
			parents.add(i);
			ranks.add(0);
		}
		this.parents = new PersistentArrayList<Integer>(parents);
		this.ranks = new PersistentArrayList<Integer>(ranks);
	}

	private PersistentPartition(PersistentArrayList<Integer> parents, PersistentArrayList<Integer> ranks) {
		this.parents = parents;
		this.ranks = ranks;
	}

	/** Find the representative for the given item */
	public int componentOf(int i) {
		final Pair<PersistentArrayList<Integer>, Integer> parentsAndComponent = auxFind(i);
		// ok to mutate parents with the compressed path. doesn't change visible behavior
		parents = parentsAndComponent.first;
		return parentsAndComponent.second;
	}

	/** Helper function for componentOf */
	private Pair<PersistentArrayList<Integer>, Integer> auxFind(int i) {
		// walk upwards until you find the canonical Component (x is canonical iff x.parent == x)
		// uses a recursive path compressing algorithm
		final Integer parent = parents.get(i);
		if (parent == i) {
			return Pair.of(parents, i);
		} else {
			final Pair<PersistentArrayList<Integer>, Integer> parentsAndComponent = auxFind(parent);
			final Integer head = parentsAndComponent.second;
			return Pair.of(parents.set(i, head), head);
		}
	}

	/** Returns a new partition in which the given components are merged */
	public PersistentPartition merge(int a, int b) {
		final int aHead = componentOf(a);
		final int bHead = componentOf(b);
		if(aHead == bHead) return this;
		// add the shorter tree underneath the taller tree
		final int aRank = ranks.get(aHead);
		final int bRank = ranks.get(bHead);
		if (aRank > bRank) {
			return new PersistentPartition(parents.set(bHead, aHead), ranks);
		} else if (bRank > aRank) {
			return new PersistentPartition(parents.set(aHead, bHead), ranks);
		} else {
			// whoops, the tree got taller
			return new PersistentPartition(parents.set(bHead, aHead), ranks.set(aHead, aRank + 1));
		}
	}

	/** Determines whether the two items are in the same component or not */
	public boolean sameComponent(int a, int b) {
		return componentOf(a) == componentOf(b);
	}
}
