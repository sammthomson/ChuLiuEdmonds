package edu.cmu.cs.lti.ark.cle;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.sun.istack.internal.Nullable;

import java.util.Collections;
import java.util.List;
import java.util.Map;

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
	public class EdgeQueue {
		private final int component;
		public final Map<Integer, ExclusiveEdge> edgesBySourceNode;

		public EdgeQueue(int component) {
			this.component = component;
			this.edgesBySourceNode = Maps.newHashMap(); // Important: don't keep this sorted. need O(1) insert time
		}

		public void addEdge(Edge edge, double weight, List<Edge> replaces) {
			// only add if source is external to SCC
			if (partition.componentOf(edge.source) == component) return;
			// only keep the best edge for each source node
			if (!edgesBySourceNode.containsKey(edge.source)
					|| weight > edgesBySourceNode.get(edge.source).weight) {
				edgesBySourceNode.put(edge.source, new ExclusiveEdge(edge, replaces, weight));
			}
		}

		@Nullable public ExclusiveEdge popBestEdge() {
			final ExclusiveEdge min = Collections.min(edgesBySourceNode.values());
			if (min == null) return null;
			edgesBySourceNode.remove(min.edge.source);
			return min;
		}
	}

	EdgeQueueMap(Partition partition) {
		this.partition = partition;
		this.queueByDestination = Maps.newHashMap();
	}

	public void addEdge(Edge edge, double weight) {
		final int destination = partition.componentOf(edge.destination);
		if (!queueByDestination.containsKey(destination)) {
			queueByDestination.put(destination, new EdgeQueue(destination));
		}
		final List<Edge> replaces = Lists.newLinkedList();
		queueByDestination.get(destination).addEdge(edge, weight, replaces);
	}

	@Nullable public ExclusiveEdge popBestEdge(int component) {
		if (!queueByDestination.containsKey(component)) return null;
		return queueByDestination.get(component).popBestEdge();
	}

	public EdgeQueue merge(int component, Iterable<Pair<EdgeQueue, Weighted<Edge>>> queuesToMerge) {
		EdgeQueue result = new EdgeQueue(component);
		for (Pair<EdgeQueue, Weighted<Edge>> queueAndReplace : queuesToMerge) {
			final EdgeQueue queue = queueAndReplace.first;
			final Weighted<Edge> replace = queueAndReplace.second;
			for (ExclusiveEdge wEdgeAndExcluded : queue.edgesBySourceNode.values()) {
				final List<Edge> replaces = wEdgeAndExcluded.excluded;
				replaces.add(replace.val);
				result.addEdge(wEdgeAndExcluded.edge, wEdgeAndExcluded.weight - replace.weight, replaces);
			}
		}
		queueByDestination.put(component, result);
		return result;
	}
}
