package edu.cmu.cs.ark.cle;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.collect.*;
import com.google.common.primitives.Doubles;

import java.util.*;

import static com.google.common.collect.ImmutableList.copyOf;
import static com.google.common.collect.Iterables.concat;
import static edu.cmu.cs.ark.cle.Weighted.weighted;
import static java.util.Collections.singleton;

/**
 * Chu-Liu-Edmonds' algorithm for finding a maximum branching in a complete, directed graph in O(n^2) time.
 * This implementation is based on Tarjan's "Finding Optimum Branchings" paper.
 * http://cw.felk.cvut.cz/lib/exe/fetch.php/courses/a4m33pal/cviceni/tarjan-finding-optimum-branchings.pdf
 *
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmonds {
	private static Function<Edge, Integer> getDestination = new Function<Edge, Integer>() {
		@Override public Integer apply(Edge input) {
			return input.destination;
		}
	};

	/** Represents the subgraph that gets iteratively built up in the CLE algorithm. */
	private static class Subgraph {
		// Partition representing the strongly connected components (SCCs).
		private Partition stronglyConnected;
		// Partition representing the weakly connected components (WCCs).
		private Partition weaklyConnected;
		// An invariant of the CLE algorithm is that each SCC always has at most one incoming edge.
		// You can think of these edges as implicitly defining a graph with SCCs as nodes.
		private final Map<Integer, Weighted<Edge>> incomingEdgeByScc;
		// History of edges we've added, and for each, a list of edges it would exclude.
		// More recently added edges get priority over less recently added edges when reconstructing the final tree.
		private final LinkedList<ExclusiveEdge> edgesAndWhatTheyExclude;
		// a priority queue of incoming edges for each SCC.
		private final EdgeQueueMap unseenIncomingEdges;
		// running sum of weights.
		// edge weights are adjusted as we go to take into account the fact that we have an extra edge in each cycle
		private double score;

		public Subgraph(double[][] graph, Integer root) {
			this(graph, root, ImmutableList.<Edge>of(), ImmutableList.<Edge>of());
		}

		public Subgraph(double[][] graph, Integer root, List<Edge> required, List<Edge> banned) {
			final int numNodes = graph.length;
			stronglyConnected = new Partition(numNodes);
			weaklyConnected = new Partition(numNodes);
			incomingEdgeByScc = Maps.newHashMap();
			edgesAndWhatTheyExclude = Lists.newLinkedList();
			unseenIncomingEdges = getEdgesByDestination(graph, root, required, banned);
			score = 0.0;
		}

		/** Groups edges by their destination component. O(n^2) */
		private EdgeQueueMap getEdgesByDestination(double[][] graph, Integer root, List<Edge> required, List<Edge> banned) {
			final ImmutableListMultimap<Integer, Edge> requiredByDestination = Multimaps.index(required, getDestination);
			final EdgeQueueMap incomingEdges = new EdgeQueueMap(stronglyConnected);
			for (int destinationNode = 0; destinationNode < graph.length; destinationNode++) {
				if(destinationNode != root) { // Throw out incoming edges for the root node.
					final ImmutableList<Edge> requiredEdgesForDest = requiredByDestination.get(destinationNode);
					final Optional<Integer> requiredDest = requiredEdgesForDest.isEmpty() ? Optional.<Integer>absent() : Optional.of(requiredEdgesForDest.get(0).source);
					for (int sourceNode = 0; sourceNode < graph.length; sourceNode++) {
						if (sourceNode == destinationNode) continue; // Skip autocycle edges
						if (requiredDest.isPresent() && sourceNode != requiredDest.get()) {
							// Skip any edge that might compete with a required edge
							continue;
						}
						if (banned.contains(new Edge(sourceNode, destinationNode))) {
							// Skip banned edges
							continue;
						}
						final double weight = graph[sourceNode][destinationNode];
						if (weight != Double.NEGATIVE_INFINITY) {
							incomingEdges.addEdge(new Edge(sourceNode, destinationNode), weight);
						}
					}
				}
			}
			return incomingEdges;
		}

		/**
		 * Given an edge that completes a cycle, merge all SCCs on that cycle into one SCC.
		 * Returns the new component.
		 */
		private int merge(Weighted<Edge> newEdge, EdgeQueueMap unseenIncomingEdges) {
			// Find edges connecting SCCs on the path from newEdge.destination to newEdge.source
			final List<Weighted<Edge>> cycle = getCycle(newEdge);
			// build up list of queues that need to be merged, with the edge they would exclude
			final List<Pair<EdgeQueueMap.EdgeQueue, Weighted<Edge>>> queuesToMerge = Lists.newLinkedList();
			for (Weighted<Edge> currentEdge : cycle) {
				final int destination = stronglyConnected.componentOf(currentEdge.val.destination);
				final EdgeQueueMap.EdgeQueue queue =
						unseenIncomingEdges.queueByDestination.get(destination);
				// if we choose an edge in queue, we'll have to throw out currentEdge at the end
				// (each SCC can have only one incoming edge).
				queuesToMerge.add(Pair.of(queue, currentEdge));
				unseenIncomingEdges.queueByDestination.remove(destination);
			}
			// Merge all SCCs on the cycle into one
			for (Weighted<Edge> e : cycle) {
				stronglyConnected.merge(e.val.source, e.val.destination);
			}
			int component = stronglyConnected.componentOf(newEdge.val.destination);
			// merge the queues and put the merged queue back into our map under the new component
			unseenIncomingEdges.merge(component, queuesToMerge);
			// keep our implicit graph of SCCs up to date:
			// we just created a cycle, so all in-edges have sources inside the new component
			// i.e. there is no edge with source outside component, and destination inside component
			incomingEdgeByScc.remove(component);
			return component;
		}

		/** Gets the cycle of edges between SCCs that newEdge creates */
		private List<Weighted<Edge>> getCycle(Weighted<Edge> newEdge) {
			final List<Weighted<Edge>> cycle = Lists.newLinkedList();
			// circle around backward in the implicit graph until you get back to where you started
			Weighted<Edge> edge = newEdge;
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
		public Optional<Integer> addEdge(ExclusiveEdge wEdgeAndExcludes) {
			final Edge edge = wEdgeAndExcludes.edge;
			final double weight = wEdgeAndExcludes.weight;
			final Weighted<Edge> wEdge = weighted(edge, weight);
			score += weight;
			final int destinationScc = stronglyConnected.componentOf(edge.destination);
			edgesAndWhatTheyExclude.addFirst(wEdgeAndExcludes);
			incomingEdgeByScc.put(destinationScc, wEdge);
			if (!weaklyConnected.sameComponent(edge.source, edge.destination)) {
				// Edge connects two different WCCs. Including it won't create a new cycle
				weaklyConnected.merge(edge.source, edge.destination);
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
			final Set<Edge> excluded = Sets.newHashSet();
			// start with the most recent
			while (!edgesAndWhatTheyExclude.isEmpty()) {
				final ExclusiveEdge edgeAndWhatItExcludes = edgesAndWhatTheyExclude.pollFirst();
				final Edge edge = edgeAndWhatItExcludes.edge;
				if(!excluded.contains(edge)) {
					excluded.addAll(edgeAndWhatItExcludes.excluded);
					parents.put(edge.destination, edge.source);
				}
			}
			return weighted(parents, score);
		}

		public Optional<ExclusiveEdge> popBestEdge(int component) {
			return popBestEdge(component, Maps.<Integer, Integer>newHashMap());
		}

		public Optional<ExclusiveEdge> popBestEdge(int component, Map<Integer, Integer> best) {
			return unseenIncomingEdges.popBestEdge(component, best);
		}
	}

	/**
	 * Find an optimal branching of the given graph, rooted in the given node.
	 * This is the main entry point for the algorithm.
	 */
	public static Weighted<Map<Integer,Integer>> getMaxSpanningTree(double[][] graph, int root) {
		return getMaxSpanningTree(graph, root, ImmutableList.<Edge>of(), ImmutableList.<Edge>of());
	}

	public static Weighted<Map<Integer,Integer>>
			getMaxSpanningTree(double[][] graph, int root, List<Edge> required, List<Edge> banned) {
		final int numNodes = graph.length;
		// result
		final Subgraph subgraph = new Subgraph(graph, root, required, banned);

		// In the beginning, subgraph has no edges, so no SCC has in-edges.
		final Queue<Integer> componentsWithNoInEdges = Lists.newLinkedList();
		for(int i = 0; i < numNodes; i++) componentsWithNoInEdges.add(i);

		// Work our way through all componentsWithNoInEdges, in no particular order
		while (!componentsWithNoInEdges.isEmpty()) {
			final int component = componentsWithNoInEdges.poll();
			// find maximum edge entering 'component' from the outside.
			final Optional<ExclusiveEdge> oMaxInEdge = subgraph.popBestEdge(component);
			if (!oMaxInEdge.isPresent()) continue; // No in-edges left to consider for this component. Done with it!
			final ExclusiveEdge maxInEdge = oMaxInEdge.get();
			// add the new edge to subgraph, merging SCCs if necessary
			final Optional<Integer> newComponent = subgraph.addEdge(maxInEdge);
			if (newComponent.isPresent()) {
				// addEdge created a cycle, which means the new cycle doesn't have any incoming edges
				componentsWithNoInEdges.add(newComponent.get());
			}
		}
		// Once no component has incoming edges left to consider, it's time to recover the optimal branching.
		return subgraph.getParentsMap();
	}

	/** Corresponds to the NEXT function in Camerini et al. 1980 */
	public static Optional<Pair<Edge, Double>> next(double[][] graph,
													int root,
													List<Edge> required,
													List<Edge> banned,
													Weighted<Map<Integer,Integer>> best) {
		final int numNodes = graph.length;
		// result
		final Subgraph subgraph = new Subgraph(graph, root, required, banned);

		// In the beginning, subgraph has no edges, so no SCC has in-edges.
		final Queue<Integer> componentsWithNoInEdges = Lists.newLinkedList();
		for(int i = 0; i < numNodes; i++) componentsWithNoInEdges.add(i);

		double bestDifference = Double.POSITIVE_INFINITY;
		Optional<ExclusiveEdge> bestEdgeToKickOut = Optional.absent();

		// Work our way through all componentsWithNoInEdges, in no particular order
		while (!componentsWithNoInEdges.isEmpty()) {
			final int component = componentsWithNoInEdges.poll();
			// find maximum edge entering 'component' from the outside.
			// break ties in favor of edges in best
			final Optional<ExclusiveEdge> oMaxInEdge = subgraph.popBestEdge(component, best.val);
			if (!oMaxInEdge.isPresent()) continue; // No in-edges left to consider for this component. Done with it!
			final ExclusiveEdge maxInEdge = oMaxInEdge.get();
			if (best.val.get(maxInEdge.edge.destination) == maxInEdge.edge.source && !required.contains(maxInEdge.edge)) {

				// TODO: left off here
				final Optional<ExclusiveEdge> oAlternativeEdge =
						seek(maxInEdge, best.val, subgraph.unseenIncomingEdges.queueByDestination.get(component));
				if (oAlternativeEdge.isPresent()) {
					final ExclusiveEdge alternativeEdge = oAlternativeEdge.get();
					final double difference = maxInEdge.weight - alternativeEdge.weight;
					if (difference < bestDifference) {
						bestDifference = difference;
						bestEdgeToKickOut = Optional.of(maxInEdge);
					}
				}
			}
			// add the new edge to subgraph, merging SCCs if necessary
			final Optional<Integer> newComponent = subgraph.addEdge(maxInEdge);
			if (newComponent.isPresent()) {
				// addEdge created a cycle, which means the new cycle doesn't have any incoming edges
				componentsWithNoInEdges.add(newComponent.get());
			}
		}
		// Once no component has incoming edges left to consider, it's time to recover the optimal branching.
		if (bestEdgeToKickOut.isPresent()) {
			return Optional.of(Pair.of(bestEdgeToKickOut.get().edge, bestDifference));
		} else {
			return Optional.absent();
		}
	}

	/** Determines whether potentialAncestor is an ancestor of node in bestArborescence */
	private static boolean isAncestor(int node, int potentialAncestor, Map<Integer, Integer> bestArborescence) {
		int currentNode = node;
		while (bestArborescence.containsKey(currentNode)) {
			currentNode = bestArborescence.get(currentNode);
			if (currentNode == potentialAncestor) return true;
		}
		return false;
	}

	/**
	 * Finds nextBestEdge, the next best alternative to maxInEdge for which the tail of
	 * maxInEdge is not an ancestor of the source of nextBestEdge in bestArborescence
	 */
	public static Optional<ExclusiveEdge>
			seek(ExclusiveEdge maxInEdge, Map<Integer, Integer> bestArborescence, EdgeQueueMap.EdgeQueue edgeQueue) {
		Optional<ExclusiveEdge> oNextBestEdge = edgeQueue.popBestEdge();
		while (oNextBestEdge.isPresent()) {
			final ExclusiveEdge nextBestEdge = oNextBestEdge.get();
			if (!isAncestor(nextBestEdge.edge.source, maxInEdge.edge.destination, bestArborescence)) {
				edgeQueue.addEdge(nextBestEdge.edge, nextBestEdge.weight, nextBestEdge.excluded);
				return oNextBestEdge;
			} else {
				oNextBestEdge = edgeQueue.popBestEdge();
			}
		}
		return Optional.absent();
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

	private static class ElementOfPartitionOfSolutions implements Comparable<ElementOfPartitionOfSolutions> {
		final double weight;
		final Edge edgeToKickOut;
		final Weighted<Map<Integer, Integer>> currentBest;
		final List<Edge> required;
		final List<Edge> banned;

		public ElementOfPartitionOfSolutions(double weight, Edge edgeToKickOut, Weighted<Map<Integer, Integer>> currentBest, List<Edge> required, List<Edge> banned) {
			this.weight = weight;
			this.edgeToKickOut = edgeToKickOut;
			this.currentBest = currentBest;
			this.required = required;
			this.banned = banned;
		}

		@Override
		public int compareTo(ElementOfPartitionOfSolutions other) {
			return Doubles.compare(weight, other.weight);
		}
	}

	public static List<Weighted<Map<Integer, Integer>>> getKBestSpanningTrees(double[][] weights, int root, int k) {
		final List<Edge> initialRequired = ImmutableList.of();
		final List<Edge> initialBanned = ImmutableList.of();
		final PriorityQueue<ElementOfPartitionOfSolutions> queue =
				new PriorityQueue<ElementOfPartitionOfSolutions>(3 * k, Collections.reverseOrder());  // want to pop max
		final List<Weighted<Map<Integer, Integer>>> results = Lists.newArrayList();
		// 1-best
		final Weighted<Map<Integer, Integer>> best = getMaxSpanningTree(weights, root);
		results.add(best);
		Optional<Pair<Edge, Double>> oEdgeToKickOutAndDiff = next(weights, root, initialRequired, initialBanned, best);
		if (oEdgeToKickOutAndDiff.isPresent()) {
			Pair<Edge, Double> edgeToKickOutAndDiff = oEdgeToKickOutAndDiff.get();
			final Edge edge = edgeToKickOutAndDiff.first;
			Double diff = edgeToKickOutAndDiff.second;
			queue.add(new ElementOfPartitionOfSolutions(
					best.weight - diff,
					edge,
					best,
					initialRequired,
					initialBanned
			));
		}
		for (int j = 1; j < k && !queue.isEmpty(); j++) {
			final ElementOfPartitionOfSolutions element = queue.poll();
			Edge edgeToKickOut = element.edgeToKickOut;
			final List<Edge> newRequired = copyOf(concat(element.required, singleton(edgeToKickOut)));
			final ImmutableList<Edge> newBanned = copyOf(concat(element.banned, singleton(edgeToKickOut)));
			final Weighted<Map<Integer, Integer>> jthBest = getMaxSpanningTree(weights, root, element.required, newBanned);
			results.add(jthBest);
			oEdgeToKickOutAndDiff = next(weights, root, newRequired, element.banned, element.currentBest);
			if (oEdgeToKickOutAndDiff.isPresent()) {
				final Pair<Edge, Double> edgeToKickOutAndDiff = oEdgeToKickOutAndDiff.get();
				final double diff = edgeToKickOutAndDiff.second;
				edgeToKickOut = edgeToKickOutAndDiff.first;
				queue.add(new ElementOfPartitionOfSolutions(
						element.currentBest.weight - diff,
						edgeToKickOut,
						element.currentBest,
						newRequired,
						element.banned
				));
			}
			oEdgeToKickOutAndDiff = next(weights, root, element.required, newBanned, jthBest);
			if (oEdgeToKickOutAndDiff.isPresent()) {
				final Pair<Edge, Double> edgeToKickOutAndDiff = oEdgeToKickOutAndDiff.get();
				final double diff = edgeToKickOutAndDiff.second;
				edgeToKickOut = edgeToKickOutAndDiff.first;
				queue.add(new ElementOfPartitionOfSolutions(
						element.weight - diff,
						edgeToKickOut,
						jthBest,
						element.required,
						newBanned
				));
			}
		}
		return results;
	}
}