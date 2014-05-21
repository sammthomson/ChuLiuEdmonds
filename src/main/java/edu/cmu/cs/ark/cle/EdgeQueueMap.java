package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import java.util.*;

/**
 * A set of queues optimized for keeping track of unprocessed Edges entering each component during the CLE algorithm.
 *
 * addEdge: O(inverseAckermann(n))
 * 	   (basically O(1), and exactly O(1) at the beginning when all components are singletons)
 * popBestEdge: O(n)
 * merge: O(n)
 * in number of nodes
 *
 * @author sthomson@cs.cmu.edu
 */
class EdgeQueueMap {
	Partition partition;
	public final Map<Integer, EdgeQueue> queueByDestination;

	/**
	 * A queue of edges entering one component.
	 * The trick is to keep track of only the best edge for each source node.
	 * This means the number of edges in the queue is bounded by the number of nodes.
	 */
	public static class EdgeQueue {
		private final int component;
		public final PriorityQueue<ExclusiveEdge> edges;
		private final Partition partition;

		public EdgeQueue(int component, Partition partition) {
			this.component = component;
			this.edges = new PriorityQueue<ExclusiveEdge>(11, Collections.reverseOrder()); // Important: don't keep this sorted. need O(1) insert time
			this.partition = partition;
		}

		public void addEdge(ExclusiveEdge exclusiveEdge) {
			// only add if source is external to SCC
			if (partition.componentOf(exclusiveEdge.edge.source) == component) return;
			// only keep the best edge for each source node
			edges.add(exclusiveEdge);
		}

		public Optional<ExclusiveEdge> popBestEdge() {
			return popBestEdge(ImmutableMap.<Integer, Integer>of());  // TODO: inefficient
		}

		/** Always breaks ties in favor of edges in bestArborescence */
		public Optional<ExclusiveEdge> popBestEdge(Map<Integer, Integer> bestArborescence) {
			final List<ExclusiveEdge> maxInEdges = maxWithTies(edges);
			if (maxInEdges.isEmpty()) return Optional.absent();
			Optional<ExclusiveEdge> maxInEdge = Optional.absent();
			for (ExclusiveEdge ee : maxInEdges) {
				final Edge e = ee.edge;
				final int dest = e.destination;
				if (bestArborescence.containsKey(dest) && bestArborescence.get(dest) == e.source) {
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

	EdgeQueueMap(Partition partition) {
		this.partition = partition;
		this.queueByDestination = Maps.newHashMap();
	}

	public void addEdge(Edge edge, double weight) {
		final int destination = partition.componentOf(edge.destination);
		if (!queueByDestination.containsKey(destination)) {
			queueByDestination.put(destination, new EdgeQueue(destination, partition));
		}
		final List<Edge> replaces = Lists.newLinkedList();
		queueByDestination.get(destination).addEdge(new ExclusiveEdge(edge, replaces, weight));
	}

	/** Always breaks ties in favor of edges in best */
	public Optional<ExclusiveEdge> popBestEdge(int component, Map<Integer, Integer> best) {
		if (!queueByDestination.containsKey(component)) return Optional.absent();
		return queueByDestination.get(component).popBestEdge(best);
	}

	public EdgeQueue merge(int component, Iterable<Pair<EdgeQueue, Weighted<Edge>>> queuesToMerge) {
		EdgeQueue result = new EdgeQueue(component, partition);
		for (Pair<EdgeQueue, Weighted<Edge>> queueAndReplace : queuesToMerge) {
			final EdgeQueue queue = queueAndReplace.first;
			final Weighted<Edge> replace = queueAndReplace.second;
			for (ExclusiveEdge wEdgeAndExcluded : queue.edges) {
				final List<Edge> replaces = wEdgeAndExcluded.excluded;
				replaces.add(replace.val);
				result.addEdge(new ExclusiveEdge(wEdgeAndExcluded.edge, replaces, wEdgeAndExcluded.weight - replace.weight));
			}
		}
		queueByDestination.put(component, result);
		return result;
	}
}
