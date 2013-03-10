package edu.cmu.cs.lti.ark.mst.graph;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import java.util.Collection;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Represents a partition of a collection into disjoint sets (components).
 */
public class Partition<T> {
	/** Map from element to the set that contains it. */
	private final ImmutableMap<T, Component> components;

	/**
	 * Represents one of the disjoint sets that make up a partition.
	 * (These are the traditional Union-find data structure, and can be used without
	 * the enclosing Partition instance if desired.)
	 */
	public static class Component {
		private Component parent;
		private int rank;

		private Component() {
			parent = this;
			rank = 0;
		}

		/** Returns the canonical Component equivalent to this. Compresses the path as it goes.  */
		private Component find() {
			// walk upwards until you find the canonical Component (x is canonical iff x.parent == x)
			// uses a single-pass path tree compressing algorithm
			Component current = this;
			Component last = this;
			while (current.parent != current) {
				last.parent = current.parent; //initially a no-op
				last = current;
				current = current.parent;
			}
			this.rank = 0;
			return current;
		}

		/** Merges other and this */
		public Component mergeWith(Component other) {
			checkNotNull(other);
			Component me = find();
			Component o = other.find();
			// add the shorter tree underneath the taller tree
			if (me == o) {
				return me; // no-op
			} else if (me.rank < o.rank) {
				me.parent = o;
				return o;
			} else if (me.rank > o.rank) {
				o.parent = me;
				return me;
			} else {
				// whoops, the tree got taller
				o.parent = me;
				me.rank++;
				return me;
			}
		}

		@Override public boolean equals(Object obj) {
			return (obj instanceof Partition.Component) && (find() == ((Partition.Component) obj).find());
		}

		// dangerous to have a mutable hashCode. be careful!
		@Override public int hashCode() {
			return System.identityHashCode(find());
		}
	}

	/** Create one component for each unique value of partitionFunction. */
	public <V> Partition(Collection<T> elements, Function<T, V> partitionFunction) {
		final Multimap<V, T> byPartition = Multimaps.index(elements, partitionFunction);
		ImmutableMap.Builder<T, Component> componentsBuilder = ImmutableMap.builder();
		for (V uniqueValue : byPartition.keySet()) {
			final Component component = new Component();
			for (T element : byPartition.get(uniqueValue)) {
				componentsBuilder.put(element, component);
			}
		}
		components = componentsBuilder.build();
	}

	/** Static constructor. Create one component for each element. */
	public static <T> Partition<T> newSingletonsPartition(Collection<T> elements) {
		return new Partition<T>(elements, Functions.<T>identity());
	}

	/** Get the set that contains element */
	public Component componentOf(T element) {
		return components.get(element);
	}

	public boolean sameComponent(T a, T b) {
		return componentOf(a).equals(componentOf(b));
	}

	/** Get all components */
	public Collection<Component> getAllComponents() {
		return components.values();
	}
}