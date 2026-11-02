package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import edu.cmu.cs.ark.cle.util.Pair;
import edu.cmu.cs.ark.cle.ds.Partition;
import edu.cmu.cs.ark.cle.graph.Edge;
import edu.cmu.cs.ark.cle.graph.WeightedGraph;
import edu.cmu.cs.ark.cle.util.Weighted;

import java.util.*;

import static com.google.common.base.Predicates.and;
import static com.google.common.base.Predicates.not;
import static edu.cmu.cs.ark.cle.EdgeQueueMap.EdgeQueue;
import static edu.cmu.cs.ark.cle.util.Weighted.weighted;

/**
 * Chu-Liu-Edmonds' algorithm for finding a maximum spanning arborescence in a directed graph.
 * Uses Tarjan's cycle contractions with the growth-path processing order of Gabow et al. (1986).
 * See https://doi.org/10.1007/BF02579168 and https://arxiv.org/abs/2208.02590, Sections 2.3 and 3.
 *
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmonds {
	/** Represents the subgraph that gets iteratively built up in the CLE algorithm. */
	static class PartialSolution<V> {
		// Partition representing the strongly connected components (SCCs).
		private final Partition<V> stronglyConnected;
		private final Iterator<V> startingNodes;
		private final Set<V> visited = Sets.newHashSet();
		private final Set<V> onPath = Sets.newHashSet();
		// Chosen edges in predecessor order; a cycle is a suffix of this path.
		private final Deque<Weighted<Edge<V>>> path = new ArrayDeque<>();
		private Optional<V> pendingComponent = Optional.absent();
		// History of edges we've added, and for each, a list of edges it would exclude.
		// More recently added edges get priority over less recently added edges when reconstructing the final tree.
		private final Deque<ExclusiveEdge<V>> edgesAndWhatTheyExclude = new ArrayDeque<>();
		// a priority queue of incoming edges for each SCC that we haven't chosen an incoming edge for yet.
		final EdgeQueueMap<V> unseenIncomingEdges;
		// running sum of weights.
		// edge weights are adjusted as we go to take into account the fact that we have an extra edge in each cycle
		private double score;

		private PartialSolution(Partition<V> stronglyConnected,
								EdgeQueueMap<V> unseenIncomingEdges,
								Iterator<V> startingNodes) {
			this.stronglyConnected = stronglyConnected;
			this.unseenIncomingEdges = unseenIncomingEdges;
			this.startingNodes = startingNodes;
		}

		public static <T> PartialSolution<T> initialize(WeightedGraph<T> graph) {
			final Partition<T> stronglyConnected = Partition.singletons(graph.getNodes());
			// group edges by their destination component
			final EdgeQueueMap<T> incomingEdges = new EdgeQueueMap<T>(stronglyConnected);
			for (T destinationNode : graph.getNodes()) {
				for (Weighted<Edge<T>> inEdge : graph.getIncomingEdges(destinationNode)) {
					if (inEdge.weight != Double.NEGATIVE_INFINITY) {
						incomingEdges.addEdge(inEdge);
					}
				}
			}
			return new PartialSolution<>(stronglyConnected, incomingEdges, graph.getNodes().iterator());
		}

		/** Continue at the last edge's source, or start a new path when the previous one is finished. */
		public Optional<V> nextComponent() {
			if (pendingComponent.isPresent()) {
				final Optional<V> next = pendingComponent;
				pendingComponent = Optional.absent();
				return next;
			}
			while (startingNodes.hasNext()) {
				final V component = stronglyConnected.componentOf(startingNodes.next());
				if (visited.add(component)) {
					path.clear();
					onPath.clear();
					onPath.add(component);
					return Optional.of(component);
				}
			}
			return Optional.absent();
		}

		/**
		 * Merge all SCCs on the given cycle into one SCC.
		 * Returns the new component.
		 */
		private V merge(List<Weighted<Edge<V>>> cycle) {
			// build up list of queues that need to be merged, with the edge they would exclude
			final List<Pair<EdgeQueue<V>, Weighted<Edge<V>>>> queuesToMerge = Lists.newArrayList();
			for (Weighted<Edge<V>> currentEdge : cycle) {
				final V destination = stronglyConnected.componentOf(currentEdge.val.destination);
				final EdgeQueue<V> queue = unseenIncomingEdges.queueByDestination.remove(destination);
				// if we choose an edge in `queue`, we'll have to throw out `currentEdge` at the end
				// (each SCC can have only one incoming edge).
				queuesToMerge.add(Pair.of(queue, currentEdge));
			}
			// Merge all SCCs on the cycle into one
			for (Weighted<Edge<V>> e : cycle) {
				stronglyConnected.merge(e.val.source, e.val.destination);
			}
			final V component = stronglyConnected.componentOf(cycle.get(0).val.destination);
			// merge the queues and put the merged queue back into our map under the new component
			unseenIncomingEdges.merge(component, queuesToMerge);
			return component;
		}

		private List<Weighted<Edge<V>>> removeCycleSuffix(V repeatedComponent) {
			final List<Weighted<Edge<V>>> cycle = Lists.newArrayList();
			V destination;
			do {
				final Weighted<Edge<V>> edge = path.removeLast();
				destination = stronglyConnected.componentOf(edge.val.destination);
				onPath.remove(destination);
				cycle.add(edge);
			} while (!destination.equals(repeatedComponent));
			return cycle;
		}

		/** Add the edge and continue the growth path, contracting its suffix if the source repeats. */
		public void addEdge(ExclusiveEdge<V> wEdgeAndExcludes) {
			final Edge<V> edge = wEdgeAndExcludes.edge;
			final Weighted<Edge<V>> wEdge = weighted(edge, wEdgeAndExcludes.weight);
			score += wEdgeAndExcludes.weight;
			edgesAndWhatTheyExclude.addFirst(wEdgeAndExcludes);
			path.addLast(wEdge);
			final V source = stronglyConnected.componentOf(edge.source);
			if (onPath.contains(source)) {
				final V component = merge(removeCycleSuffix(source));
				onPath.add(component);
				pendingComponent = Optional.of(component);
			} else if (visited.add(source)) {
				onPath.add(source);
				pendingComponent = Optional.of(source);
			}
		}

		/**
		 * Recovers the optimal arborescence.
		 *
		 * Each SCC can only have 1 edge entering it: the edge that we added most recently.
		 * So we work backwards, adding edges unless they conflict with edges we've already added.
		 * Runtime is O(n^2) in the worst case.
		 */
		private Weighted<Arborescence<V>> recoverBestArborescence() {
			final ImmutableMap.Builder<V, V> parents = ImmutableMap.builder();
			final Set<Edge<V>> excluded = Sets.newHashSet();
			// start with the most recent
			while (!edgesAndWhatTheyExclude.isEmpty()) {
				final ExclusiveEdge<V> edgeAndWhatItExcludes = edgesAndWhatTheyExclude.removeFirst();
				final Edge<V> edge = edgeAndWhatItExcludes.edge;
				if (!excluded.contains(edge)) {
					excluded.addAll(edgeAndWhatItExcludes.excluded);
					parents.put(edge.destination, edge.source);
				}
			}
			return weighted(Arborescence.of(parents.build()), score);
		}

		public Optional<ExclusiveEdge<V>> popBestEdge(V component) {
			return popBestEdge(component, Arborescence.<V>empty());
		}

		/** Always breaks ties in favor of edges in `best` */
		public Optional<ExclusiveEdge<V>> popBestEdge(V component, Arborescence<V> best) {
			return unseenIncomingEdges.popBestEdge(component, best);
		}
	}

	/**
	 * Find an optimal arborescence of the given graph `graph`, rooted in the given node `root`.
	 */
	public static <V> Weighted<Arborescence<V>> getMaxArborescence(WeightedGraph<V> graph, V root) {
		// remove all edges incoming to `root`. resulting arborescence is then forced to be rooted at `root`.
		return getMaxArborescence(graph.filterEdges(not(Edge.hasDestination(root))));
	}

	static <V> Weighted<Arborescence<V>> getMaxArborescence(WeightedGraph<V> graph,
															Set<Edge<V>> required,
															Set<Edge<V>> banned) {
		return getMaxArborescence(graph.filterEdges(and(not(Edge.competesWith(required)), not(Edge.isIn(banned)))));
	}

	/**
	 * Find an optimal arborescence of the given graph.
	 */
	public static <V> Weighted<Arborescence<V>> getMaxArborescence(WeightedGraph<V> graph) {
		final PartialSolution<V> partialSolution =
				PartialSolution.initialize(graph.filterEdges(not(Edge.<V>isAutoCycle())));
		for (Optional<V> next = partialSolution.nextComponent(); next.isPresent();
			 next = partialSolution.nextComponent()) {
			final V component = next.get();
			// find maximum edge entering 'component' from outside 'component'.
			final Optional<ExclusiveEdge<V>> oMaxInEdge = partialSolution.popBestEdge(component);
			if (!oMaxInEdge.isPresent()) continue; // No in-edges left to consider for this component. Done with it!
			partialSolution.addEdge(oMaxInEdge.get());
		}
		// Once no component has incoming edges left to consider, it's time to recover the optimal branching.
		return partialSolution.recoverBestArborescence();
	}
}
