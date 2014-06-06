package edu.cmu.cs.ark.cle;

import com.google.common.base.Preconditions;
import com.google.common.collect.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import static com.google.common.collect.DiscreteDomain.integers;
import static com.google.common.collect.Range.closedOpen;

/**
 * @author sthomson@cs.cmu.edu
 */
public class DenseWeightedGraph<V> implements WeightedGraph<V> {
	final private ArrayList<V> nodes;
	final private Map<V, Integer> indexOf;
	final private double[][] weights;

	private DenseWeightedGraph(ArrayList<V> nodes, Map<V, Integer> indexOf, double[][] weights) {
		this.nodes = nodes;
		this.indexOf = indexOf;
		this.weights = weights;
	}

	public static <V> DenseWeightedGraph<V> from(Iterable<V> nodes, double[][] weights) {
		final ArrayList<V> nodeList = Lists.newArrayList(nodes);
		Preconditions.checkArgument(nodeList.size() == weights.length);
		final Map<V, Integer> indexOf = Maps.newHashMap();
		for (int i = 0; i < nodeList.size(); i++) {
			indexOf.put(nodeList.get(i), i);
		}
		return new DenseWeightedGraph<V>(nodeList, indexOf, weights);
	}

	public static DenseWeightedGraph<Integer> from(double[][] weights) {
		final ContiguousSet<Integer> nodes = ContiguousSet.create(closedOpen(0, weights.length), integers());
		return DenseWeightedGraph.from(nodes, weights);
	}

	@Override
	public Collection<V> getNodes() {
		return nodes;
	}

	@Override
	public double getWeightOf(V source, V dest) {
		if (!indexOf.containsKey(source) || !indexOf.containsKey(dest)) return Double.NEGATIVE_INFINITY;
		return weights[indexOf.get(source)][indexOf.get(dest)];
	}
}
