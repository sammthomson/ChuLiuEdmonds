package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

/**
 * @author sthomson@cs.cmu.edu
 */
public class SparseWeightedGraph<V> extends WeightedGraph<V> {
	final private Set<V> nodes;
	final private Map<V, Map<V, Double>> outgoingEdges;
	final private Map<V, Map<V, Weighted<Edge<V>>>> incomingEdges;

	private SparseWeightedGraph(Set<V> nodes,
								Map<V, Map<V, Double>> outgoingEdges,
								Map<V, Map<V, Weighted<Edge<V>>>> incomingEdges) {
		this.nodes = nodes;
		this.outgoingEdges = outgoingEdges;
		this.incomingEdges = incomingEdges;
	}

	public static <T> SparseWeightedGraph<T> from(Iterable<T> nodes, Iterable<Weighted<Edge<T>>> edges) {
		final Map<T, Map<T, Double>> outgoingEdges = Maps.newHashMap();
		final Map<T, Map<T, Weighted<Edge<T>>>> incomingEdges = Maps.newHashMap();
		for (Weighted<Edge<T>> edge : edges) {
			if (!outgoingEdges.containsKey(edge.val.source)) {
				outgoingEdges.put(edge.val.source, Maps.<T, Double>newHashMap());
			}
			outgoingEdges.get(edge.val.source).put(edge.val.destination, edge.weight);
			if (!incomingEdges.containsKey(edge.val.destination)) {
				incomingEdges.put(edge.val.destination, Maps.<T, Weighted<Edge<T>>>newHashMap());
			}
			incomingEdges.get(edge.val.destination).put(edge.val.source, edge);
		}
		return new SparseWeightedGraph<T>(ImmutableSet.copyOf(nodes), outgoingEdges, incomingEdges);
	}

	public static <T> SparseWeightedGraph<T> from(Iterable<Weighted<Edge<T>>> edges) {
		final Set<T> nodes = Sets.newHashSet();
		for (Weighted<Edge<T>> edge : edges) {
			nodes.add(edge.val.source);
			nodes.add(edge.val.destination);
		}
		return SparseWeightedGraph.from(nodes, edges);
	}

	@Override
	public Collection<V> getNodes() {
		return nodes;
	}

	@Override
	public double getWeightOf(V source, V dest) {
		if (!outgoingEdges.containsKey(source)) return Double.NEGATIVE_INFINITY;
		final Map<V, Double> outEdges = outgoingEdges.get(source);
		if (!outEdges.containsKey(dest)) return Double.NEGATIVE_INFINITY;
		return outEdges.get(dest);
	}

	@Override
	public Collection<Weighted<Edge<V>>> getIncomingEdges(V destinationNode) {
		if (!incomingEdges.containsKey(destinationNode)) return ImmutableSet.of();
		return incomingEdges.get(destinationNode).values();
	}
}
