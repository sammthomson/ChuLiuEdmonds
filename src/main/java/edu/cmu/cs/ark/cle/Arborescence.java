package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableMap;

import java.util.Map;

/**
 * @author sthomson@cs.cmu.edu
 */
public class Arborescence<V> {
	/**
	 * In an arborescence, each node (other than the root) has exactly one parent. This is the map
	 * from each node to its parent.
	 */
	final Map<V, V> parents;

	private Arborescence(Map<V, V> parents) {
		this.parents = parents;
	}

	public static <T> Arborescence<T> of(Map<T, T> parents) {
		return new Arborescence<T>(parents);
	}

	public static <T> Arborescence<T> empty() {
		return Arborescence.of(ImmutableMap.<T, T>of());
	}
}
