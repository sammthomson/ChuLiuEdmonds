package edu.cmu.cs.lti.ark.cle;

import com.google.common.collect.Maps;
import com.sun.istack.internal.Nullable;

import java.util.Collections;
import java.util.Map;

import static edu.cmu.cs.lti.ark.cle.Weighted.weighted;

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
	PersistentPartition partition;
	public final Map<Integer, EdgeQueue> queueByDestination;

	/**
	 * A queue of edges entering one component.
	 * The trick is to keep track of only the best edge for each source node.
	 * This means the number of edges in the queue is bounded by the number of nodes.
	 */
	public class EdgeQueue {
		private final int component;
		public final Map<Integer, Weighted<Edge<Integer>>> edgesBySourceNode;

		public EdgeQueue(int component) {
			this.component = component;
			this.edgesBySourceNode = Maps.newHashMap(); // Important: don't keep this sorted. need O(1) insert time
		}

		public void addEdge(Weighted<Edge<Integer>> wEdge) {
			final Edge<Integer> edge = wEdge.val;
			// only add if source is external to SCC
			if (partition.componentOf(edge.source) == component) return;
			// only keep the best edge for each source node
			if (!edgesBySourceNode.containsKey(edge.source)
					|| wEdge.compareTo(edgesBySourceNode.get(edge.source)) < 0) {
				edgesBySourceNode.put(edge.source, wEdge);
			}
		}

		@Nullable public Weighted<Edge<Integer>> popBestEdge() {
			final Weighted<Edge<Integer>> min = Collections.min(edgesBySourceNode.values());
			if (min == null) return null;
			edgesBySourceNode.remove(min.val.source);
			return min;
		}
	}

	EdgeQueueMap(PersistentPartition partition) {
		this.partition = partition;
		this.queueByDestination = Maps.newHashMap();
	}

	public void addEdge(Weighted<Edge<Integer>> wEdge) {
		final Integer destination = partition.componentOf(wEdge.val.destination);
		if (!queueByDestination.containsKey(destination)) {
			queueByDestination.put(destination, new EdgeQueue(destination));
		}
		queueByDestination.get(destination).addEdge(wEdge);
	}

	@Nullable public Weighted<Edge<Integer>> popBestEdge(int component) {
		if (!queueByDestination.containsKey(component)) return null;
		return queueByDestination.get(component).popBestEdge();
	}

	public EdgeQueue merge(int component, Iterable<Weighted<EdgeQueue>> queuesToMerge) {
		EdgeQueue result = new EdgeQueue(component);
		for (Weighted<EdgeQueue> q : queuesToMerge) {
			for (Weighted<Edge<Integer>> wEdge : q.val.edgesBySourceNode.values()) {
				result.addEdge(doOffset(wEdge, q.weight));
			}
		}
		queueByDestination.put(component, result);
		return result;
	}

	/** Add offset to the weight of wEdge */
	private Weighted<Edge<Integer>> doOffset(Weighted<Edge<Integer>> wEdge, double offset) {
		return weighted(new Edge<Integer>(wEdge.val.source, wEdge.val.destination), wEdge.weight + offset);
	}
}
