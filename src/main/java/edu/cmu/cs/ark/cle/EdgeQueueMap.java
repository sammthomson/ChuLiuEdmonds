package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;


class EdgeQueueMap<V> {
	final Partition<V> partition;
	public final Map<V, EdgeQueue<V>> queueByDestination;

	public static class EdgeQueue<V> {
		private final V component;
		public final PriorityQueue<ExclusiveEdge<V>> edges;
		private final Partition<V> partition;

		private EdgeQueue(V component, Partition<V> partition) {
			this.component = component;
			this.edges = new PriorityQueue<ExclusiveEdge<V>>(11, Collections.reverseOrder());
			this.partition = partition;
		}

		public static <T> EdgeQueue<T> create(T component, Partition<T> partition) {
			return new EdgeQueue<T>(component, partition);
		}

		public void addEdge(ExclusiveEdge<V> exclusiveEdge) {
			// only add if source is external to SCC
			if (partition.componentOf(exclusiveEdge.edge.source) == component) return;
			// only keep the best edge for each source node
			edges.add(exclusiveEdge);
		}

		public Optional<ExclusiveEdge<V>> popBestEdge() {
			return popBestEdge(Arborescence.<V>empty());
		}

		/** Always breaks ties in favor of edges in bestArborescence */
		public Optional<ExclusiveEdge<V>> popBestEdge(Arborescence<V> bestArborescence) {
			final List<ExclusiveEdge<V>> maxInEdges = maxWithTies(edges);
			if (maxInEdges.isEmpty()) return Optional.absent();
			Optional<ExclusiveEdge<V>> maxInEdge = Optional.absent();
			for (ExclusiveEdge<V> ee : maxInEdges) {
				final Edge<V> e = ee.edge;
				final V dest = e.destination;
				if (bestArborescence.parents.containsKey(dest) && bestArborescence.parents.get(dest) == e.source) {
					maxInEdge = Optional.of(ee);
					break;
				}
			}
			if (!maxInEdge.isPresent()) {
				maxInEdge = Optional.of(maxInEdges.get(0));
			}
			edges.remove(maxInEdge.get());
			return maxInEdge;
		}
	}

	public static <T extends Object & Comparable<? super T>> List<T> maxWithTies(PriorityQueue<T> queue) {
		List<T> candidates = Lists.newArrayList();
		if (!queue.isEmpty()) {
			candidates.add(queue.poll());
		}
		while (!queue.isEmpty()) {
			T next = queue.poll();
			final int cmp = queue.comparator().compare(candidates.get(0), next);
			if (cmp == 0) {
				candidates.add(next);
			} else {
				queue.add(next);
				break;
			}
		}
		for (T candidate : candidates) {
			queue.add(candidate);
		}
		return candidates;
	}

	EdgeQueueMap(Partition<V> partition) {
		this.partition = partition;
		this.queueByDestination = Maps.newHashMap();
	}

	public void addEdge(Edge<V> edge, double weight) {
		final V destination = partition.componentOf(edge.destination);
		if (!queueByDestination.containsKey(destination)) {
			queueByDestination.put(destination, EdgeQueue.create(destination, partition));
		}
		final List<Edge<V>> replaces = Lists.newLinkedList();
		queueByDestination.get(destination).addEdge(new ExclusiveEdge<V>(edge, replaces, weight));
	}

	/** Always breaks ties in favor of edges in best */
	public Optional<ExclusiveEdge<V>> popBestEdge(V component, Arborescence<V> best) {
		if (!queueByDestination.containsKey(component)) return Optional.absent();
		return queueByDestination.get(component).popBestEdge(best);
	}

	public EdgeQueue merge(V component, Iterable<Pair<EdgeQueue<V>, Weighted<Edge<V>>>> queuesToMerge) {
		EdgeQueue<V> result = EdgeQueue.create(component, partition);
		for (Pair<EdgeQueue<V>, Weighted<Edge<V>>> queueAndReplace : queuesToMerge) {
			final EdgeQueue<V> queue = queueAndReplace.first;
			final Weighted<Edge<V>> replace = queueAndReplace.second;
			for (ExclusiveEdge<V> wEdgeAndExcluded : queue.edges) {
				final List<Edge<V>> replaces = wEdgeAndExcluded.excluded;
				replaces.add(replace.val);
				result.addEdge(new ExclusiveEdge<V>(wEdgeAndExcluded.edge, replaces, wEdgeAndExcluded.weight - replace.weight));
			}
		}
		queueByDestination.put(component, result);
		return result;
	}
}
