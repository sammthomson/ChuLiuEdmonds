package edu.cmu.cs.lti.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.*;
import com.sun.istack.internal.Nullable;

import java.util.*;

import static edu.cmu.cs.lti.ark.cle.Partition.Component;
import static edu.cmu.cs.lti.ark.cle.Weighted.weighted;

/**
 * Chu-Liu-Edmonds' algorithm for finding a maximum branching in a complete, directed graph.
 * This implementation is based on Tarjan's "Finding Optimum Branchings" paper.
 * http://cw.felk.cvut.cz/lib/exe/fetch.php/courses/a4m33pal/cviceni/tarjan-finding-optimum-branchings.pdf
 *
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmonds {
	/**
	 * Represents the subgraph that gets iteratively built up in the CLE algorithm.
	 * @param <T> The node type
	 */
	private static class Subgraph<T> {
		private static final int EXPECTED_BRANCH_FACTOR = 5; // not that important
		// Number of nodes in the graph
		private final int numNodes;
		// Map from node to the edges that emanate from it
		private final Multimap<T, Edge<T>> edgesBySource;
		// Partition representing the strongly connected components (SCCs).
		private final Partition<T> stronglyConnected;
		// Partition representing the weakly connected components (WCCs).
		private final Partition<T> weaklyConnected;
		// An invariant of the CLE algorithm is that each SCC always has at most one incoming edge.
		private final Map<Component, Weighted<Edge<T>>> incomingEdgeByScc;
		// running sum of weights. takes into account the fact that we have an extra edge for each cycle
		private double score;

		public Subgraph(Collection<T> nodes) {
			numNodes = nodes.size();
			edgesBySource = HashMultimap.create(numNodes, EXPECTED_BRANCH_FACTOR);
			stronglyConnected = Partition.newSingletonsPartition(nodes);
			weaklyConnected = Partition.newSingletonsPartition(nodes);
			incomingEdgeByScc = Maps.newHashMap();
			score = 0.0;
		}

		/** Gets the maximum weight incomingEdge whose source is actually outside component */
		@Nullable
		private Weighted<Edge<T>> getMaxIncomingEdge(Component component, Queue<Weighted<Edge<T>>> incomingEdges) {
			Weighted<Edge<T>> edge = incomingEdges.poll();
			while(edge != null && stronglyConnected.componentOf(edge.val.source).equals(component)) {
				edge = incomingEdges.poll(); // internal edge, keep trying
			}
			return edge;
		}

		/** Given an edge that completes a cycle, merge all SCCs on that cycle into one SCC. */
		private Component merge(Weighted<Edge<T>> newEdge, Map<Component, PriorityQueue<Weighted<Edge<T>>>> unseenIncomingEdges) {
			// Find edges connecting SCCs on the path from maxInEdge.destination to maxInEdge.source
			final List<Weighted<Edge<T>>> cycle = getCycle(newEdge);
			// get the minimum edge weight on the cycle (edges are naturally ordered max-first)
			//final double minEdgeWeight = Collections.min(cycle, Ordering.natural().reverse()).getWeight(); // this is what the paper says, but doesn't make sense to me
			// Increment edge weights in each queue and then add them to the merged SCC queue.
			final PriorityQueue<Weighted<Edge<T>>> mergedQueue = Queues.newPriorityQueue();
			for (Weighted<Edge<T>> currentEdge : cycle) {
				//final double inc = minEdgeWeight - currentEdge.getWeight(); // this is what the paper says, but doesn't make sense and doesn't work
				final double inc = -currentEdge.weight;
				final Component destination = stronglyConnected.componentOf(currentEdge.val.destination);
				for (Weighted<Edge<T>> wEdge : unseenIncomingEdges.get(destination)) {
					Edge<T> e = wEdge.val;
					mergedQueue.add(weighted(new Edge<T>(e.source, e.destination), wEdge.weight + inc));
				}
				unseenIncomingEdges.remove(destination);
			}
			// Merge all SCCs on the cycle into one
			Component component = stronglyConnected.componentOf(newEdge.val.destination);
			for (Weighted<Edge<T>> e : cycle) {
				component = component.mergeWith(stronglyConnected.componentOf(e.val.source));
			}
			// put the merged queue back into our map
			unseenIncomingEdges.put(component, mergedQueue);
			// we just created a cycle, so all in-edges have sources inside the component now
			// i.e. there is no edge with source outside component, and destination inside component
			incomingEdgeByScc.remove(component);
			return component;
		}

		/** Gets the cycle of edges between SCCs that newEdge creates */
		private List<Weighted<Edge<T>>> getCycle(Weighted<Edge<T>> newEdge) {
			final List<Weighted<Edge<T>>> cycle = Lists.newArrayList();
			// circle around until you get back to where you started
			Weighted<Edge<T>> edge = newEdge;
			cycle.add(edge);
			while (!stronglyConnected.sameComponent(edge.val.source, newEdge.val.destination)) {
				edge = incomingEdgeByScc.get(stronglyConnected.componentOf(edge.val.source));
				cycle.add(edge);
			}
			return cycle;
		}

		/**
		 * Adds the given edge to this subgraph, merging SCCs if necessary
		 * @return the new SCC, if adding edge created a cycle
		 */
		public Optional<Component> addEdge(Weighted<Edge<T>> wEdge, Map<Component, PriorityQueue<Weighted<Edge<T>>>> unseenIncomingEdges) {
			Edge<T> edge = wEdge.val;
			edgesBySource.put(edge.source, edge);
			score += wEdge.weight;
			final Component sourceWcc = weaklyConnected.componentOf(edge.source);
			final Component destinationWcc = weaklyConnected.componentOf(edge.destination);
			if (!sourceWcc.equals(destinationWcc)) {
				// Edge connects two different WCCs. Including it won't create a new cycle
				sourceWcc.mergeWith(destinationWcc);
				incomingEdgeByScc.put(stronglyConnected.componentOf(edge.destination), wEdge);
				return Optional.absent();
			} else {
				// Edge is contained within one WCC. Including it will create a new cycle.
				return Optional.of(merge(wEdge, unseenIncomingEdges));
			}
		}

		/**
		 * Gets the optimal spanning tree, encoded as a map from each node to its parent.
		 *
		 * Since what we have is a bunch of cycles each with <= 1 incoming edge,
		 * our strategy is to follow edges forward from the root, but throw out the
		 * last edge in each cycle.
		 */
		private Weighted<Map<T, T>> getParentsMap(T root) {
			final Map<T, T> parents = Maps.newHashMap();
			// keep track of what we've seen, so we can throw out the last edge in each cycle
			final Set<T> visited = Sets.newHashSetWithExpectedSize(numNodes);
			final LinkedList<T> toDo = Lists.newLinkedList();
			toDo.add(root);
			while(!toDo.isEmpty()) {
				T node = toDo.pop();
				visited.add(node);
				for (Edge<T> outEdge : edgesBySource.get(node)) {
					if (!visited.contains(outEdge.destination)) {
						parents.put(outEdge.destination, node);
						toDo.add(outEdge.destination);
					}
				}
			}
			return weighted(parents, score);
		}
	}

	/**
	 * Find an optimal branching of the given graph, rooted in the given node.
	 * This is the main entry point for the algorithm.
	 */
	public static Weighted<Map<Integer, Integer>> getMaxSpanningTree(double[][] originalGraph, int root) {
		final int numNodes = originalGraph.length;
		List<Integer> nodes = Lists.newArrayListWithExpectedSize(numNodes);
		for(int i = 0; i < numNodes; i++) nodes.add(i);
		// result
		final Subgraph<Integer> subgraph = new Subgraph<Integer>(nodes);

		// a priority queue of incoming edges for each SCC.
		final Map<Component, PriorityQueue<Weighted<Edge<Integer>>>> unseenIncomingEdges =
				getEdgesByDestination(originalGraph, root, subgraph.stronglyConnected);
		// In the beginning, subgraph has no edges, so no SCC has in-edges.
		final LinkedList<Component> componentsWithNoInEdges =
				Lists.newLinkedList(subgraph.stronglyConnected.getAllComponents());

		// Work our way through all components, in no particular order
		while (!componentsWithNoInEdges.isEmpty()) {
			final Component component = componentsWithNoInEdges.pop();

			// find maximum edge entering 'component' from the outside.
			final Weighted<Edge<Integer>> maxInEdge =
					subgraph.getMaxIncomingEdge(component, unseenIncomingEdges.get(component));
			if (maxInEdge == null) continue; // No in-edges left to consider for this component. Done with it!
			// add the new edge to subgraph, merging SCCs if necessary
			final Optional<Component> newScc = subgraph.addEdge(maxInEdge, unseenIncomingEdges);
			if (newScc.isPresent()) {
				// addEdge created a cycle, which means the new cycle doesn't have any incoming edges
				componentsWithNoInEdges.add(newScc.get());
			}
		}
		// Once no component has incoming edges left to consider, we can
		// recover the optimal branching.
		return subgraph.getParentsMap(root);
	}

	/**
	 * Groups edges by their destination component
	 *
	 * TODO: this is O(m^2 log(m)), need to use a specialized PriorityQueue to get down to O(m^2)
	 */
	private static Map<Component, PriorityQueue<Weighted<Edge<Integer>>>>
			getEdgesByDestination(double[][] originalGraph, Integer root, Partition<Integer> scc) {
		final int numNodes = originalGraph.length;
		final Map<Component, PriorityQueue<Weighted<Edge<Integer>>>> incomingEdges = Maps.newHashMapWithExpectedSize(numNodes);
		for (int destinationNode = 0; destinationNode < numNodes; destinationNode++) {
			// Create a priority queue of incoming edges for each SCC.
			final PriorityQueue<Weighted<Edge<Integer>>> sccPriorityQueue = Queues.newPriorityQueue();
			if(destinationNode != root) { // Throw out incoming edges for the root node.
				for (int sourceNode = 0; sourceNode < numNodes; sourceNode++) {
					if (sourceNode == destinationNode) continue; // Skip autocycle edges
					final double weight = originalGraph[sourceNode][destinationNode];
					if (weight != Double.NEGATIVE_INFINITY)
						sccPriorityQueue.add(weighted(new Edge<Integer>(sourceNode, destinationNode), weight));
				}
			}
			incomingEdges.put(scc.componentOf(destinationNode), sccPriorityQueue);
		}
		return incomingEdges;
	}
}