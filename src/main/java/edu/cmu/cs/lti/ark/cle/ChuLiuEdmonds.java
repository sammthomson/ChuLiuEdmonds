package edu.cmu.cs.lti.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.*;

import java.util.*;

import static edu.cmu.cs.lti.ark.cle.Partition.Component;
import static edu.cmu.cs.lti.ark.cle.Weighted.weighted;

/**
 * Chu-Liu-Edmonds' algorithm for finding a maximum branching in a complete, directed graph in O(n^2) time.
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
		// You can think of these edges as defining a graph with SCCs as nodes.
		private final Map<Component, Weighted<Edge<T>>> incomingEdgeByScc;
		// running sum of weights.
		// edge weights are adjusted as we go to take into account the fact that we have an extra edge in each cycle
		private double score;

		public Subgraph(Collection<T> nodes) {
			numNodes = nodes.size();
			edgesBySource = HashMultimap.create(numNodes, EXPECTED_BRANCH_FACTOR);
			stronglyConnected = Partition.newSingletonsPartition(nodes);
			weaklyConnected = Partition.newSingletonsPartition(nodes);
			incomingEdgeByScc = Maps.newHashMap();
			score = 0.0;
		}

		/** Given an edge that completes a cycle, merge all SCCs on that cycle into one SCC. */
		private Component merge(Weighted<Edge<T>> newEdge, EdgeQueueMap<T> unseenIncomingEdges) {
			// Find edges connecting SCCs on the path from newEdge.destination to newEdge.source
			final List<Weighted<Edge<T>>> cycle = getCycle(newEdge);
			// build up list of queues that need to be merged, with their respective weight offsets
			final List<Weighted<EdgeQueueMap.EdgeQueue<T>>> queuesToMerge = Lists.newLinkedList();
			for (Weighted<Edge<T>> currentEdge : cycle) {
				final Component destination = stronglyConnected.componentOf(currentEdge.val.destination);
				final EdgeQueueMap.EdgeQueue<T> queue =
						unseenIncomingEdges.queueByDestination.get(destination);
				// if we choose an edge in queue, we'll have to throw out currentEdge at the end
				// (each node can have only one parent).
				// offset the weight of every edge in queue to reflect that
				queuesToMerge.add(weighted(queue, -currentEdge.weight));
				unseenIncomingEdges.queueByDestination.remove(destination);
			}
			// Merge all SCCs on the cycle into one
			Component component = stronglyConnected.componentOf(newEdge.val.destination);
			for (Weighted<Edge<T>> e : cycle) {
				component = component.mergeWith(stronglyConnected.componentOf(e.val.source));
			}
			// merge the queues and put the merged queue back into our map under the new component
			unseenIncomingEdges.merge(component, queuesToMerge);
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
		public Optional<Component> addEdge(Weighted<Edge<T>> wEdge, EdgeQueueMap<T> unseenIncomingEdges) {
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
			final Queue<T> toDo = Lists.newLinkedList();
			toDo.add(root);
			while(!toDo.isEmpty()) {
				T node = toDo.poll();
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
	public static Weighted<Map<Integer,Integer>> getMaxSpanningTree(double[][] graph, int root) {
		final int numNodes = graph.length;
		List<Integer> nodes = Lists.newArrayListWithExpectedSize(numNodes);
		for(int i = 0; i < numNodes; i++) nodes.add(i);
		// result
		final Subgraph<Integer> subgraph = new Subgraph<Integer>(nodes);

		// a priority queue of incoming edges for each SCC.
		final EdgeQueueMap<Integer> unseenIncomingEdges =
				getEdgesByDestination(graph, root, subgraph);
		// In the beginning, subgraph has no edges, so no SCC has in-edges.
		final Queue<Component> componentsWithNoInEdges =
				Lists.newLinkedList(subgraph.stronglyConnected.getAllComponents());

		// Work our way through all components, in no particular order
		while (!componentsWithNoInEdges.isEmpty()) {
			final Component component = componentsWithNoInEdges.poll();

			// find maximum edge entering 'component' from the outside.
			final Weighted<Edge<Integer>> maxInEdge = unseenIncomingEdges.popBestEdge(component);
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
	 * Find the k best branchings of the given graph, rooted in the given node.
	 * Induce diversity by penalizing results that share edges with previous results.
	 *
	 * @param originalGraph the graph to find branchings for
	 * @param root which node the branchings must be rooted on
	 * @param k number of best branchings to return
	 * @param alpha the factor by which to penalize repeated edges
	 * @return a list of the k best branchings, along with their scores
	 */
	public static List<Weighted<Map<Integer,Integer>>> getDiverseKBestSpanningTrees(double[][] originalGraph,
																					int root,
																					int k,
																					double alpha) {
		final double[][] graph = copyGraph(originalGraph); // we're about to mutate this
		final List<Weighted<Map<Integer, Integer>>> results = Lists.newArrayListWithExpectedSize(k);
		for (int i = 0; i < k; i++) {
			final Weighted<Map<Integer, Integer>> maxSpanningTree = getMaxSpanningTree(graph, root);
			results.add(maxSpanningTree);
			// penalize edges for appearing in a previous solution
			for (int to : maxSpanningTree.val.keySet()) {
				int from = maxSpanningTree.val.get(to);
				graph[from][to] -= alpha;
			}
		}
		return results;
	}

	private static double[][] copyGraph(double[][] originalGraph) {
		final double[][] result = new double[originalGraph.length][];
		for (int i = 0; i < originalGraph.length; i++) result[i] = originalGraph[i].clone();
		return result;
	}

	/** Groups edges by their destination component. O(n^2) */
	private static EdgeQueueMap<Integer> getEdgesByDestination(double[][] graph,
															   Integer root,
															   Subgraph<Integer> subgraph) {
		final EdgeQueueMap<Integer> incomingEdges = new EdgeQueueMap<Integer>(subgraph.stronglyConnected);
		for (int destinationNode = 0; destinationNode < graph.length; destinationNode++) {
			if(destinationNode != root) { // Throw out incoming edges for the root node.
				for (int sourceNode = 0; sourceNode < graph.length; sourceNode++) {
					if (sourceNode == destinationNode) continue; // Skip autocycle edges
					final double weight = graph[sourceNode][destinationNode];
					if (weight != Double.NEGATIVE_INFINITY)
						incomingEdges.addEdge(weighted(new Edge<Integer>(sourceNode, destinationNode), weight));
				}
			}
		}
		return incomingEdges;
	}
}