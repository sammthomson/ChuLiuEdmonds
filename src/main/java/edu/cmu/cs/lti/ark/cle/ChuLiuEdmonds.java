package edu.cmu.cs.lti.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.*;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

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
		// A stack of our history of SCCs over time. Needed to reconstruct the tree at the end.
		// PersistentPartitions share structure, so doesn't take up much space.
		private final LinkedList<PersistentPartition> sccHistory;
		// Partition representing the weakly connected components (WCCs).
		private PersistentPartition weaklyConnected;
		// An invariant of the CLE algorithm is that each SCC always has at most one incoming edge.
		// You can think of these edges as defining a graph with SCCs as nodes.
		private final Map<Integer, Weighted<Edge<Integer>>> incomingEdgeByCurrentScc;
		// Edges we've added, arranged in rows by history, and columns by the SCC they enter
		// Also used when reconstructing the final tree
		private final Table<PersistentPartition, Integer, Weighted<Edge<Integer>>> incomingEdgeByHistoricScc;
		// running sum of weights.
		// edge weights are adjusted as we go to take into account the fact that we have an extra edge in each cycle
		private double score;

		public Subgraph(int numNodes) {
			stronglyConnected = new PersistentPartition(numNodes);
			sccHistory = Lists.newLinkedList();
			weaklyConnected = new PersistentPartition(numNodes);
			incomingEdgeByHistoricScc = HashBasedTable.create();
			incomingEdgeByCurrentScc = Maps.newHashMap();
			score = 0.0;
		}

		/** Given an edge that completes a cycle, merge all SCCs on that cycle into one SCC. */
		private PersistentPartition merge(Weighted<Edge<Integer>> newEdge, EdgeQueueMap unseenIncomingEdges) {
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
			sccHistory.addFirst(stronglyConnected);
			for (Weighted<Edge<Integer>> e : cycle) {
				stronglyConnected = stronglyConnected.merge(e.val.source, e.val.destination);
			}
			int component = stronglyConnected.componentOf(newEdge.val.destination);
			// merge the queues and put the merged queue back into our map under the new component
			unseenIncomingEdges.partition = stronglyConnected;
			unseenIncomingEdges.merge(component, queuesToMerge);
			// we just created a cycle, so all in-edges have sources inside the component now
			// i.e. there is no edge with source outside component, and destination inside component
			incomingEdgeByCurrentScc.remove(component);
			return stronglyConnected;
		}

		/** Gets the cycle of edges between SCCs that newEdge creates */
		private List<Weighted<Edge<Integer>>> getCycle(Weighted<Edge<Integer>> newEdge) {
			final List<Weighted<Edge<Integer>>> cycle = Lists.newArrayList();
			// circle around until you get back to where you started
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
		 * @return the new Partition, if adding edge created a cycle
		 */
		public Optional<PersistentPartition> addEdge(Weighted<Edge<Integer>> wEdge, EdgeQueueMap unseenIncomingEdges) {
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
		 * Since what we have is a bunch of cycles each with <= 1 incoming edge,
		 * our strategy is to follow edges forward from the root, but throw out the
		 * last edge in each cycle.
		 */
		private Weighted<Map<Integer, Integer>> getParentsMap() {
			final Map<Integer, Integer> parents = Maps.newHashMap();
			sccHistory.addFirst(stronglyConnected);
			// unpeel history, layer by layer
			while (!sccHistory.isEmpty()) {
				final PersistentPartition historicPartition = sccHistory.pop();
				final Map<Integer, Weighted<Edge<Integer>>> incomingEdgeByScc =
						incomingEdgeByHistoricScc.row(historicPartition);
				for(Weighted<Edge<Integer>> edge : incomingEdgeByScc.values()) {
					parents.put(edge.val.destination, edge.val.source);
					// for any SCC that this edge enters, it can be the only one
					for (PersistentPartition histScc : sccHistory) {
						// TODO: going through the PersistentPartition's history back and forth, over and over like this
						// might be inefficient with the PersistentArrayList implementation we're using
						incomingEdgeByHistoricScc.remove(histScc, histScc.componentOf(edge.val.destination));
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
		final Subgraph subgraph = new Subgraph(numNodes);

		// a priority queue of incoming edges for each SCC.
		final EdgeQueueMap unseenIncomingEdges =
				getEdgesByDestination(graph, root, subgraph);
		// In the beginning, subgraph has no edges, so no SCC has in-edges.
		final Queue<Integer> componentsWithNoInEdges = Lists.newLinkedList(nodes);

		// Work our way through all components, in no particular order
		while (!componentsWithNoInEdges.isEmpty()) {
			final int component = componentsWithNoInEdges.poll();

			// find maximum edge entering 'component' from the outside.
			final Weighted<Edge<Integer>> maxInEdge = unseenIncomingEdges.popBestEdge(component);
			if (maxInEdge == null) continue; // No in-edges left to consider for this component. Done with it!
			// add the new edge to subgraph, merging SCCs if necessary
			final Optional<PersistentPartition> newPartition = subgraph.addEdge(maxInEdge, unseenIncomingEdges);
			if (newPartition.isPresent()) {
				// addEdge created a cycle, which means the new cycle doesn't have any incoming edges
				componentsWithNoInEdges.add(newPartition.get().componentOf(component));
			}
		}
		// Once no component has incoming edges left to consider, we can
		// recover the optimal branching.
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