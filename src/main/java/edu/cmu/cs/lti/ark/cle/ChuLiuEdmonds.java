package edu.cmu.cs.lti.ark.cle;

import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.*;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import static com.google.common.collect.Iterables.any;
import static edu.cmu.cs.lti.ark.cle.Weighted.weighted;

/**
 * Chu-Liu-Edmonds' algorithm for finding a maximum branching in a complete, directed graph in O(n^2) time.
 * This implementation is based on Tarjan's "Finding Optimum Branchings" paper.
 * http://cw.felk.cvut.cz/lib/exe/fetch.php/courses/a4m33pal/cviceni/tarjan-finding-optimum-branchings.pdf
 *
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmonds {
	/** Represents the subgraph that gets iteratively built up in the CLE algorithm. */
	private static class Subgraph {
		// Partition representing the strongly connected components (SCCs).
		private PersistentPartition stronglyConnected;
		// Partition representing the weakly connected components (WCCs).
		private PersistentPartition weaklyConnected;
		// An invariant of the CLE algorithm is that each SCC always has at most one incoming edge.
		// You can think of these edges as implicitly defining a graph with SCCs as nodes.
		private final Map<Integer, Weighted<Edge<Integer>>> incomingEdgeByCurrentScc;
		// History of edges we've added. Needed to reconstruct the final tree.
		// rows are a time-slice in history, columns are the SCC they enter, values are the edge
		// PersistentPartitions share structure, and we use a sparse table implementation, so doesn't take up much space.
		private final Table<PersistentPartition, Integer, Weighted<Edge<Integer>>> incomingEdgeByHistoricScc;
		// running sum of weights.
		// edge weights are adjusted as we go to take into account the fact that we have an extra edge in each cycle
		private double score;

		public Subgraph(int numNodes) {
			stronglyConnected = new PersistentPartition(numNodes);
			weaklyConnected = new PersistentPartition(numNodes);
			incomingEdgeByCurrentScc = Maps.newHashMap();
			// keep rows sorted by natural ordering (partitions with fewest components (i.e. most recent) come first)
			incomingEdgeByHistoricScc = TreeBasedTable.create();
			score = 0.0;
		}

		/**
		 * Given an edge that completes a cycle, merge all SCCs on that cycle into one SCC.
		 * Returns the new component.
		 */
		private int merge(Weighted<Edge<Integer>> newEdge, EdgeQueueMap unseenIncomingEdges) {
			// Find edges connecting SCCs on the path from newEdge.destination to newEdge.source
			final List<Weighted<Edge<Integer>>> cycle = getCycle(newEdge);
			// build up list of queues that need to be merged, with their respective weight offsets
			final List<Weighted<EdgeQueueMap.EdgeQueue>> queuesToMerge = Lists.newLinkedList();
			for (Weighted<Edge<Integer>> currentEdge : cycle) {
				final int destination = stronglyConnected.componentOf(currentEdge.val.destination);
				final EdgeQueueMap.EdgeQueue queue =
						unseenIncomingEdges.queueByDestination.get(destination);
				// if we choose an edge in queue, we'll have to throw out currentEdge at the end
				// (each SCC can have only one incoming edge).
				// so offset the weight of every edge in queue to reflect that
				queuesToMerge.add(weighted(queue, -currentEdge.weight));
				unseenIncomingEdges.queueByDestination.remove(destination);
			}
			// Merge all SCCs on the cycle into one
			for (Weighted<Edge<Integer>> e : cycle) {
				stronglyConnected = stronglyConnected.merge(e.val.source, e.val.destination);
			}
			int component = stronglyConnected.componentOf(newEdge.val.destination);
			// merge the queues and put the merged queue back into our map under the new component
			unseenIncomingEdges.partition = stronglyConnected;
			unseenIncomingEdges.merge(component, queuesToMerge);
			// keep our implicit graph of SCCs up to date:
			// we just created a cycle, so all in-edges have sources inside the new component
			// i.e. there is no edge with source outside component, and destination inside component
			incomingEdgeByCurrentScc.remove(component);
			return stronglyConnected.componentOf(component);
		}

		/** Gets the cycle of edges between SCCs that newEdge creates */
		private List<Weighted<Edge<Integer>>> getCycle(Weighted<Edge<Integer>> newEdge) {
			final List<Weighted<Edge<Integer>>> cycle = Lists.newArrayList();
			// circle around backward in the implicit graph until you get back to where you started
			Weighted<Edge<Integer>> edge = newEdge;
			cycle.add(edge);
			while (!stronglyConnected.sameComponent(edge.val.source, newEdge.val.destination)) {
				edge = incomingEdgeByCurrentScc.get(stronglyConnected.componentOf(edge.val.source));
				cycle.add(edge);
			}
			return cycle;
		}

		/**
		 * Adds the given edge to this subgraph, merging SCCs if necessary
		 * @return the new SCC, if adding edge created a cycle
		 */
		public Optional<Integer> addEdge(Weighted<Edge<Integer>> wEdge, EdgeQueueMap unseenIncomingEdges) {
			score += wEdge.weight;
			final Edge<Integer> edge = wEdge.val;
			final int destinationScc = stronglyConnected.componentOf(edge.destination);
			incomingEdgeByHistoricScc.put(stronglyConnected, destinationScc, wEdge);
			incomingEdgeByCurrentScc.put(destinationScc, wEdge);
			if (!weaklyConnected.sameComponent(edge.source, edge.destination)) {
				// Edge connects two different WCCs. Including it won't create a new cycle
				weaklyConnected = weaklyConnected.merge(edge.source, edge.destination);
				return Optional.absent();
			} else {
				// Edge is contained within one WCC. Including it will create a new cycle.
				return Optional.of(merge(wEdge, unseenIncomingEdges));
			}
		}

		/**
		 * Gets the optimal spanning tree, encoded as a map from each node to its parent.
		 *
		 * Each SCC can only have 1 edge entering it: the edge that we added most recently.
		 * So we work backwards, adding edges unless they conflict with edges we've already added.
		 * O(n log n) (number of edges in subgraph * number of partitions in our history)
		 */
		private Weighted<Map<Integer, Integer>> getParentsMap() {
			final Map<Integer, Integer> parents = Maps.newHashMap();
			// unpeel history, layer by layer
			// start with the most recent (i.e. most merged / fewest components)
			LinkedList<PersistentPartition> sccHistory = Lists.newLinkedList(incomingEdgeByHistoricScc.rowKeySet());
			while (!sccHistory.isEmpty()) {
				final PersistentPartition historicPartition = sccHistory.pollFirst();
				for(final Weighted<Edge<Integer>> edge : incomingEdgeByHistoricScc.row(historicPartition).values()) {
					// make sure we don't already have an edge coming into the same SCC
					final Predicate<Integer> sameSccAsEdge = new Predicate<Integer>() {
						@Override public boolean apply(Integer destination) {
							return historicPartition.sameComponent(destination, edge.val.destination);
						}
					};
					if (!any(parents.keySet(), sameSccAsEdge)) {
						parents.put(edge.val.destination, edge.val.source);
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
		final List<Integer> nodes = Lists.newArrayListWithExpectedSize(numNodes);
		for(int i = 0; i < numNodes; i++) nodes.add(i);
		// result
		final Subgraph subgraph = new Subgraph(numNodes);

		// a priority queue of incoming edges for each SCC.
		final EdgeQueueMap unseenIncomingEdges = getEdgesByDestination(graph, root, subgraph);
		// In the beginning, subgraph has no edges, so no SCC has in-edges.
		final Queue<Integer> componentsWithNoInEdges = Lists.newLinkedList(nodes);

		// Work our way through all componentsWithNoInEdges, in no particular order
		while (!componentsWithNoInEdges.isEmpty()) {
			final int component = componentsWithNoInEdges.poll();
			// find maximum edge entering 'component' from the outside.
			final Weighted<Edge<Integer>> maxInEdge = unseenIncomingEdges.popBestEdge(component);
			if (maxInEdge == null) continue; // No in-edges left to consider for this component. Done with it!
			// add the new edge to subgraph, merging SCCs if necessary
			final Optional<Integer> newComponent = subgraph.addEdge(maxInEdge, unseenIncomingEdges);
			if (newComponent.isPresent()) {
				// addEdge created a cycle, which means the new cycle doesn't have any incoming edges
				componentsWithNoInEdges.add(newComponent.get());
			}
		}
		// Once no component has incoming edges left to consider, it's time to recover the optimal branching.
		return subgraph.getParentsMap();
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
	public static List<Weighted<Map<Integer,Integer>>>
			getDiverseKBestSpanningTrees(double[][] originalGraph, int root, int k, double alpha) {
		// make a copy; we're about to mutate this
		final double[][] graph = new double[originalGraph.length][];
		for (int i = 0; i < originalGraph.length; i++) graph[i] = originalGraph[i].clone();

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

	/** Groups edges by their destination component. O(n^2) */
	private static EdgeQueueMap getEdgesByDestination(double[][] graph, Integer root, Subgraph subgraph) {
		final EdgeQueueMap incomingEdges = new EdgeQueueMap(subgraph.stronglyConnected);
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