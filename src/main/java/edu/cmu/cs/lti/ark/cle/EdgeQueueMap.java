package edu.cmu.cs.lti.ark.cle;

import com.google.common.collect.Maps;
import com.sun.istack.internal.Nullable;

import java.util.Collections;
import java.util.Map;

import static edu.cmu.cs.lti.ark.cle.Weighted.weighted;

/**
 * A set of queues optimized for keeping track of unprocessed Edges entering each Component during the CLE algorithm.
 *
 * addEdge: O(inverseAckermann(n))
 * 	   (basically O(1), and exactly O(1) at the beginning when all components are singletons)
 * popBestEdge: O(n)
 * merge: O(n)
 * in number of nodes
 *
 * @param <T> the node type
 * @author sthomson@cs.cmu.edu
 */
class EdgeQueueMap<T> {
	private final Partition<T> stronglyConnected;
	public final Map<Partition.Component, EdgeQueue<T>> queueByDestination;

	/**
	 * A queue of edges entering one Component.
	 * The trick is to keep track of only the best edge for each source node.
	 * This means the number of edges in the queue is bounded by the number of nodes.
	 * @param <T> the node type
	 */
	public static class EdgeQueue<T> {
		private final Partition<T> stronglyConnected;
		private final Partition.Component component;
		public final Map<T, Weighted<Edge<T>>> edgesBySourceNode;

		public EdgeQueue(Partition<T> stronglyConnected, Partition.Component component) {
			this.stronglyConnected = stronglyConnected;
			this.component = component;
			this.edgesBySourceNode = Maps.newHashMap(); // Important: don't keep this sorted. need O(1) insert time
		}

		public void addEdge(Weighted<Edge<T>> wEdge) {
			final Edge<T> edge = wEdge.val;
			// only add if source is external to SCC
			if (stronglyConnected.componentOf(edge.source).equals(component)) return;
			// only keep the best edge for each source node
			if (!edgesBySourceNode.containsKey(edge.source)
					|| wEdge.compareTo(edgesBySourceNode.get(edge.source)) < 0) {
				edgesBySourceNode.put(edge.source, wEdge);
			}
		}

		@Nullable public Weighted<Edge<T>> popBestEdge() {
			final Weighted<Edge<T>> min = Collections.min(edgesBySourceNode.values());
			if (min == null) return null;
			edgesBySourceNode.remove(min.val.source);
			return min;
		}
	}

	EdgeQueueMap(Partition<T> stronglyConnected) {
		this.stronglyConnected = stronglyConnected;
		this.queueByDestination = Maps.newHashMap();
	}

	public void addEdge(Weighted<Edge<T>> wEdge) {
		final Partition.Component destination = stronglyConnected.componentOf(wEdge.val.destination);
		if (!queueByDestination.containsKey(destination)) {
			queueByDestination.put(destination, new EdgeQueue<T>(stronglyConnected, destination));
		}
		queueByDestination.get(destination).addEdge(wEdge);
	}

	@Nullable public Weighted<Edge<T>> popBestEdge(Partition.Component component) {
		if (!queueByDestination.containsKey(component)) return null;
		return queueByDestination.get(component).popBestEdge();
	}

	public EdgeQueue<T> merge(Partition.Component component, Iterable<Weighted<EdgeQueue<T>>> queuesToMerge) {
		EdgeQueue<T> result = new EdgeQueue<T>(stronglyConnected, component);
		for (Weighted<EdgeQueue<T>> q : queuesToMerge) {
			for (Weighted<Edge<T>> wEdge : q.val.edgesBySourceNode.values()) {
				result.addEdge(doOffset(wEdge, q.weight));
			}
		}
		queueByDestination.put(component, result);
		return result;
	}

	/** Add offset to the weight of wEdge */
	private Weighted<Edge<T>> doOffset(Weighted<Edge<T>> wEdge, double offset) {
		return weighted(new Edge<T>(wEdge.val.source, wEdge.val.destination), wEdge.weight + offset);
	}
}
