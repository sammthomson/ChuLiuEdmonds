package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableList;
import com.google.common.primitives.Doubles;
import edu.cmu.cs.ark.cle.graph.Edge;

import java.util.AbstractList;
import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * An edge, together with an immutable list of edges that can't be in the final answer if 'edge' is.
 *
 * @author sthomson@cs.cmu.edu
 */
public class ExclusiveEdge<V> implements Comparable<ExclusiveEdge<V>> {
	public final Edge<V> edge;
	public final List<Edge<V>> excluded;
	public final double weight;
	final Exclusions<V> exclusions;

	private ExclusiveEdge(Edge<V> edge, Exclusions<V> excluded, double weight) {
		this.edge = edge;
		this.excluded = excluded;
		this.exclusions = excluded;
		this.weight = weight;
	}

	public static <T> ExclusiveEdge<T> of(Edge<T> edge, List<Edge<T>> excluded, double weight) {
		return of(edge, new Exclusions<>(excluded), weight);
	}

	static <T> ExclusiveEdge<T> of(Edge<T> edge, Exclusions<T> excluded, double weight) {
		return new ExclusiveEdge<>(edge, excluded, weight);
	}

	public static <T> ExclusiveEdge<T> of(Edge<T> edge, double weight) {
		return ExclusiveEdge.of(edge, ImmutableList.<Edge<T>>of(), weight);
	}

	@Override public int compareTo(ExclusiveEdge<V> exclusiveEdge) {
		return Doubles.compare(weight, exclusiveEdge.weight);
	}

	/** Immutable concatenation tree: contractions share history instead of copying it into every edge. */
	static final class Exclusions<V> extends AbstractList<Edge<V>> {
		private final List<Edge<V>> leaf;
		private final Exclusions<V> left;
		private final Exclusions<V> right;
		private final int size;

		Exclusions(List<Edge<V>> edges) {
			leaf = ImmutableList.copyOf(edges);
			left = right = null;
			size = leaf.size();
		}

		private Exclusions(Exclusions<V> left, Exclusions<V> right) {
			this.left = left;
			this.right = right;
			leaf = null;
			size = Math.addExact(left.size, right.size);
		}

		Exclusions<V> concat(Exclusions<V> other) {
			if (isEmpty()) return other;
			if (other.isEmpty()) return this;
			return new Exclusions<>(this, other);
		}

		@Override public int size() {
			return size;
		}

		@Override public Edge<V> get(int index) {
			if (index < 0 || index >= size) throw new IndexOutOfBoundsException(index);
			Exclusions<V> current = this;
			while (current.leaf == null) {
				if (index < current.left.size) {
					current = current.left;
				} else {
					index -= current.left.size;
					current = current.right;
				}
			}
			return current.leaf.get(index);
		}

		@Override public Iterator<Edge<V>> iterator() {
			return new Iterator<>() {
				private final Deque<Exclusions<V>> pending = new ArrayDeque<>();
				private Iterator<Edge<V>> current = Collections.emptyIterator();
				{
					pending.add(Exclusions.this);
				}

				@Override public boolean hasNext() {
					while (!current.hasNext() && !pending.isEmpty()) {
						final Exclusions<V> next = pending.removeLast();
						if (next.leaf != null) {
							current = next.leaf.iterator();
						} else {
							pending.addLast(next.right);
							pending.addLast(next.left);
						}
					}
					return current.hasNext();
				}

				@Override public Edge<V> next() {
					if (!hasNext()) throw new NoSuchElementException();
					return current.next();
				}
			};
		}
	}
}
