package edu.cmu.cs.ark.cle;

import edu.cmu.cs.ark.cle.graph.Edge;
import edu.cmu.cs.ark.cle.graph.SparseWeightedGraph;
import edu.cmu.cs.ark.cle.util.Weighted;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static edu.cmu.cs.ark.cle.util.Weighted.weighted;

/** A deliberately contraction-heavy workload, not a general-purpose performance comparison. */
public final class ContractionBenchmark {
	static SparseWeightedGraph<Integer> graph(int size) {
		final List<Integer> nodes = new ArrayList<>();
		final List<Weighted<Edge<Integer>>> edges = new ArrayList<>();
		for (int node = 0; node <= 2 * size; node++) nodes.add(node);
		edges.add(weighted(Edge.from(nodes.get(0)).to(nodes.get(1)), -10.0 * size));
		for (int node = 2; node <= size; node++) {
			edges.add(weighted(Edge.from(nodes.get(node)).to(nodes.get(1)), -(double) node));
			edges.add(weighted(Edge.from(nodes.get(1)).to(nodes.get(node)), 0.0));
		}
		// These unused alternatives stay queued throughout the nested contractions.
		for (int node = size + 1; node < nodes.size(); node++) {
			edges.add(weighted(Edge.from(nodes.get(node)).to(nodes.get(1)), -20.0 * size - node));
			edges.add(weighted(Edge.from(nodes.get(0)).to(nodes.get(node)), 0.0));
		}
		return SparseWeightedGraph.from(nodes, edges);
	}

	public static void main(String[] args) {
		for (int warmup = 0; warmup < 5; warmup++) {
			ChuLiuEdmonds.getMaxArborescence(graph(100), 0);
		}
		for (int size : new int[] {500, 1000, 2000}) {
			final var graph = graph(size);
			final long[] times = new long[3];
			for (int run = 0; run < times.length; run++) {
				final long start = System.nanoTime();
				final var result = ChuLiuEdmonds.getMaxArborescence(graph, 0);
				times[run] = System.nanoTime() - start;
				if (result.weight != -10.0 * size || result.val.parents.size() != 2 * size) {
					throw new AssertionError("Incorrect benchmark result: " + result);
				}
			}
			Arrays.sort(times);
			System.out.printf("%d nodes, %d edges: median %.1f ms%n",
					2 * size + 1, 4 * size - 1, times[1] / 1e6);
		}
	}
}
