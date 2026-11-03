package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.Maps;
import edu.cmu.cs.ark.cle.ds.Partition;
import edu.cmu.cs.ark.cle.graph.Edge;
import edu.cmu.cs.ark.cle.util.Pair;
import edu.cmu.cs.ark.cle.util.Weighted;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static edu.cmu.cs.ark.cle.ExclusiveEdge.Exclusions;
import static java.util.Collections.singletonList;

class EdgeQueueMap<V> {
	final Partition<V> partition;
	public final Map<V, EdgeQueue<V>> queueByDestination;

	public static class EdgeQueue<V> {
		private final V component;
		private final Partition<V> partition;
		private HeapNode<V> root;

		private EdgeQueue(V component, Partition<V> partition) {
			this.component = component;
			this.partition = partition;
		}

		public static <T> EdgeQueue<T> create(T component, Partition<T> partition) {
			return new EdgeQueue<T>(component, partition);
		}

		public void addEdge(ExclusiveEdge<V> exclusiveEdge) {
			// only add if source is external to SCC
			if (partition.sameComponent(exclusiveEdge.edge.source, component)) return;
			root = meld(root, new HeapNode<>(exclusiveEdge));
		}

		public Optional<ExclusiveEdge<V>> popBestEdge() {
			return popBestEdge(Arborescence.<V>empty());
		}

		/** Always breaks ties in favor of edges in bestArborescence */
		public Optional<ExclusiveEdge<V>> popBestEdge(Arborescence<V> bestArborescence) {
			discardInternalEdges();
			if (root == null) return Optional.absent();
			final List<ExclusiveEdge<V>> candidates = new ArrayList<>();
			ExclusiveEdge<V> bestEdge = removeBest();
			while (!bestArborescence.contains(bestEdge.edge)) {
				discardInternalEdges();
				if (root == null || Double.compare(bestEdge.weight, root.weight) != 0) break;
				candidates.add(bestEdge);
				bestEdge = removeBest();
			}
			for (ExclusiveEdge<V> candidate : candidates) {
				addEdge(candidate);
			}
			return Optional.of(bestEdge);
		}

		private void discardInternalEdges() {
			while (root != null && partition.sameComponent(root.edge.source, component)) {
				root.push();
				root = meld(root.left, root.right);
			}
		}

		private ExclusiveEdge<V> removeBest() {
			final ExclusiveEdge<V> best = ExclusiveEdge.of(root.edge, root.exclusions, root.weight);
			root.push();
			root = meld(root.left, root.right);
			return best;
		}
	}

	/**
	 * Leftist max-heap; its right spine has logarithmic depth.
	 * A root stores its absolute weight, while a child stores its nonnegative gap below its parent,
	 * so adjusting a root implicitly adjusts its whole subtree. Since {@code parent - gap <= parent}
	 * under any rounding, floating-point error can never place a child above its parent.
	 */
	private static final class HeapNode<V> {
		private final Edge<V> edge;
		private Exclusions<V> exclusions;
		private double weight;
		private HeapNode<V> left;
		private HeapNode<V> right;
		private int rank = 1;
		private Exclusions<V> pendingExclusions;

		private HeapNode(ExclusiveEdge<V> edge) {
			this.edge = edge.edge;
			exclusions = edge.exclusions;
			weight = edge.weight;
		}

		/** Adjusts the weight of this root, and implicitly of its subtree. */
		private void adjust(double delta, Exclusions<V> excluded) {
			final double adjusted = weight + delta;
			if (Double.isFinite(weight) && !Double.isFinite(adjusted)) {
				throw new ArithmeticException("Adjusted edge weight overflows: " + weight + " + " + delta);
			}
			weight = adjusted;
			appendExclusions(excluded);
		}

		private void appendExclusions(Exclusions<V> excluded) {
			exclusions = exclusions.concat(excluded);
			pendingExclusions = pendingExclusions == null ? excluded : pendingExclusions.concat(excluded);
		}

		/** Converts this root's children from gaps to absolute weights, so they can become roots. */
		private void push() {
			materialize(left);
			materialize(right);
			pendingExclusions = null;
		}

		private void materialize(HeapNode<V> child) {
			if (child == null) return;
			child.weight = weight - child.weight;
			if (pendingExclusions != null) child.appendExclusions(pendingExclusions);
		}

		/** Converts a child from an absolute weight to its gap below this node. */
		private void rebase(HeapNode<V> child) {
			if (child == null) return;
			final double absolute = child.weight;
			// Equal weights, including equal infinities, have no gap (and must not produce NaN).
			final double gap = weight == absolute ? 0.0 : weight - absolute;
			if (!(gap >= 0.0) || (gap == Double.POSITIVE_INFINITY && Double.isFinite(absolute))) {
				throw new ArithmeticException(
						"Edge weights " + weight + " and " + absolute + " are too far apart to compare");
			}
			child.weight = gap;
		}
	}

	private static <V> HeapNode<V> meld(HeapNode<V> a, HeapNode<V> b) {
		if (a == null) return b;
		if (b == null) return a;
		if (Double.compare(a.weight, b.weight) < 0) {
			final HeapNode<V> swap = a;
			a = b;
			b = swap;
		}
		a.push();
		a.right = meld(a.right, b);
		if (a.left == null || a.left.rank < a.right.rank) {
			final HeapNode<V> swap = a.left;
			a.left = a.right;
			a.right = swap;
		}
		a.rank = a.right == null ? 1 : a.right.rank + 1;
		a.rebase(a.left);
		a.rebase(a.right);
		return a;
	}

	EdgeQueueMap(Partition<V> partition) {
		this.partition = partition;
		this.queueByDestination = Maps.newHashMap();
	}

	public void addEdge(Weighted<Edge<V>> edge) {
		final V destination = partition.componentOf(edge.val.destination);
		if (!queueByDestination.containsKey(destination)) {
			queueByDestination.put(destination, EdgeQueue.create(destination, partition));
		}
		queueByDestination.get(destination).addEdge(ExclusiveEdge.of(edge.val, edge.weight));
	}

	/** Always breaks ties in favor of edges in best */
	public Optional<ExclusiveEdge<V>> popBestEdge(V component, Arborescence<V> best) {
		if (!queueByDestination.containsKey(component)) return Optional.absent();
		return queueByDestination.get(component).popBestEdge(best);
	}

	/** Consumes the old queues; adjusting each old root implicitly adjusts every edge below it. */
	public EdgeQueue<V> merge(V component, Iterable<Pair<EdgeQueue<V>, Weighted<Edge<V>>>> queuesToMerge) {
		final EdgeQueue<V> result = EdgeQueue.create(component, partition);
		for (Pair<EdgeQueue<V>, Weighted<Edge<V>>> queueAndReplace : queuesToMerge) {
			final EdgeQueue<V> queue = queueAndReplace.first;
			final Weighted<Edge<V>> replace = queueAndReplace.second;
			if (queue.root == null) continue;
			queue.root.adjust(-replace.weight, new Exclusions<>(singletonList(replace.val)));
			result.root = meld(result.root, queue.root);
			queue.root = null;
		}
		queueByDestination.put(component, result);
		return result;
	}
}
