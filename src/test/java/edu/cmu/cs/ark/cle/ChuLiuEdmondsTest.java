package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import edu.cmu.cs.ark.cle.graph.DenseWeightedGraph;
import edu.cmu.cs.ark.cle.graph.Edge;
import edu.cmu.cs.ark.cle.graph.SparseWeightedGraph;
import edu.cmu.cs.ark.cle.graph.WeightedGraph;
import edu.cmu.cs.ark.cle.util.Weighted;
import org.junit.Test;

import java.util.*;

import static edu.cmu.cs.ark.cle.util.Weighted.weighted;
import static org.junit.Assert.*;

/**
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmondsTest {
	final static double DELTA = 0.001;
	final static double NINF = Double.NEGATIVE_INFINITY;
	final static WeightedGraph<Integer> graph = SparseWeightedGraph.from(ImmutableList.of(
			weighted(Edge.from(0).to(1), 5),
			weighted(Edge.from(0).to(2), 1),
			weighted(Edge.from(0).to(3), 1),
			weighted(Edge.from(1).to(2), 11),
			weighted(Edge.from(1).to(3), 4),
			weighted(Edge.from(2).to(1), 10),
			weighted(Edge.from(2).to(3), 5),
			weighted(Edge.from(3).to(1), 9),
			weighted(Edge.from(3).to(2), 8)
	));

	static <V> void assertEdgesSumToScore(WeightedGraph<V> originalEdgeWeights, Weighted<Arborescence<V>> bestTree) {
		final Map<V, V> parentsMap = bestTree.val.parents;
		double sumOfWeights = 0.0;
		for (V dest : parentsMap.keySet()) {
			final V source = parentsMap.get(dest);
			sumOfWeights += originalEdgeWeights.getWeightOf(source, dest);
		}
		assertEquals(sumOfWeights, bestTree.weight, DELTA);
	}

	static <V> void assertRootedTree(WeightedGraph<V> graph, V root, Weighted<Arborescence<V>> tree) {
		assertEquals(graph.getNodes().size() - 1, tree.val.parents.size());
		assertFalse(tree.val.parents.containsKey(root));
		for (V node : graph.getNodes()) {
			final Set<V> seen = new HashSet<>();
			V current = node;
			while (!current.equals(root)) {
				assertTrue("Cycle in " + tree, seen.add(current));
				assertTrue("Missing parent for " + current, tree.val.parents.containsKey(current));
				current = tree.val.parents.get(current);
			}
		}
		assertEdgesSumToScore(graph, tree);
	}

	// Independent oracle: enumerate parent assignments and reject those that do not reach the root.
	static List<Weighted<Arborescence<Integer>>> enumerateArborescences(WeightedGraph<Integer> graph, int root) {
		final List<Integer> destinations = new ArrayList<>(graph.getNodes());
		destinations.remove(Integer.valueOf(root));
		final List<Weighted<Arborescence<Integer>>> trees = new ArrayList<>();
		enumerateArborescences(graph, root, destinations, 0, new HashMap<Integer, Integer>(), 0.0, trees);
		return trees;
	}

	private static void enumerateArborescences(WeightedGraph<Integer> graph, int root,
											 List<Integer> destinations, int index, Map<Integer, Integer> parents,
											 double score, List<Weighted<Arborescence<Integer>>> trees) {
		if (index == destinations.size()) {
			for (int node : destinations) {
				final Set<Integer> seen = new HashSet<>();
				int current = node;
				while (current != root) {
					if (!seen.add(current)) return;
					current = parents.get(current);
				}
			}
			trees.add(weighted(Arborescence.of(ImmutableMap.copyOf(parents)), score));
			return;
		}
		final int destination = destinations.get(index);
		for (Weighted<Edge<Integer>> edge : graph.getIncomingEdges(destination)) {
			if (edge.val.source.equals(destination) || edge.weight == NINF) continue;
			parents.put(destination, edge.val.source);
			enumerateArborescences(graph, root, destinations, index + 1, parents, score + edge.weight, trees);
		}
		parents.remove(destination);
	}

	static WeightedGraph<Integer> nestedCyclesGraph() {
		return SparseWeightedGraph.from(ImmutableList.of(6, 5, 4, 3, 2, 1, 0), ImmutableList.of(
				weighted(Edge.from(0).to(1), 1),
				weighted(Edge.from(1).to(2), 10),
				weighted(Edge.from(2).to(1), 9),
				weighted(Edge.from(2).to(3), 8),
				weighted(Edge.from(3).to(2), 7),
				weighted(Edge.from(3).to(4), 6),
				weighted(Edge.from(4).to(3), 5),
				weighted(Edge.from(0).to(4), 0),
				weighted(Edge.from(4).to(5), 4),
				weighted(Edge.from(5).to(6), 3),
				weighted(Edge.from(0).to(6), -1)));
	}

	@Test
	public void testNestedCyclesWithPathLeadingIntoThem() {
		final WeightedGraph<Integer> graph = nestedCyclesGraph();
		final Weighted<Arborescence<Integer>> tree = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		assertRootedTree(graph, 0, tree);
		assertEquals(32.0, tree.weight, DELTA);
	}

	@Test
	public void testManyNestedCyclesWithUnusedIncomingEdges() {
		final int size = 128;
		final WeightedGraph<Integer> graph = ContractionBenchmark.graph(size);
		final Weighted<Arborescence<Integer>> tree = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		final ImmutableMap.Builder<Integer, Integer> parents = ImmutableMap.builder();
		parents.put(1, 0);
		for (int node = 2; node <= 2 * size; node++) {
			parents.put(node, node <= size ? 1 : 0);
		}
		assertEquals(parents.build(), tree.val.parents);
		assertEquals(-10.0 * size, tree.weight, DELTA);
		assertRootedTree(graph, 0, tree);
	}

	@Test
	public void testPositiveInfiniteWeightsCompeteWithFiniteWeights() {
		final WeightedGraph<Integer> forced = SparseWeightedGraph.from(
				ImmutableList.of(0, 1, 2), ImmutableList.of(
						weighted(Edge.from(0).to(1), Double.POSITIVE_INFINITY),
						weighted(Edge.from(2).to(1), 1.0),
						weighted(Edge.from(0).to(2), 1.0)));
		final Weighted<Arborescence<Integer>> tree = ChuLiuEdmonds.getMaxArborescence(forced, 0);
		assertEquals(ImmutableMap.of(1, 0, 2, 0), tree.val.parents);
		assertEquals(Double.POSITIVE_INFINITY, tree.weight, 0.0);

		// Contracting the cycle 1 <-> 2 subtracts +Infinity from the other edges into node 1.
		final WeightedGraph<Integer> cyclic = SparseWeightedGraph.from(
				ImmutableList.of(0, 1, 2), ImmutableList.of(
						weighted(Edge.from(2).to(1), Double.POSITIVE_INFINITY),
						weighted(Edge.from(1).to(2), 5.0),
						weighted(Edge.from(0).to(1), 1.0),
						weighted(Edge.from(0).to(2), 2.0)));
		final Weighted<Arborescence<Integer>> cyclicTree = ChuLiuEdmonds.getMaxArborescence(cyclic, 0);
		assertEquals(ImmutableMap.of(1, 2, 2, 0), cyclicTree.val.parents);
	}

	@Test
	public void testMixedScaleWeightsDoNotSelectSuboptimalTree() {
		// Each node in {2, 3, 4} must be reached through node 1, so the optimal tree takes 0 -> 1 (0.5).
		// Rounding adjusted weights to different values once picked 5 -> 1 (0.25) instead.
		final WeightedGraph<Integer> graph = SparseWeightedGraph.from(
				ImmutableList.of(0, 1, 2, 3, 4, 5), ImmutableList.of(
						weighted(Edge.from(2).to(1), 1e16),
						weighted(Edge.from(4).to(1), 1.875),
						weighted(Edge.from(3).to(1), 1.625),
						weighted(Edge.from(0).to(1), 0.5),
						weighted(Edge.from(5).to(1), 0.25),
						weighted(Edge.from(0).to(5), 0),
						weighted(Edge.from(1).to(3), 0),
						weighted(Edge.from(1).to(4), 0),
						weighted(Edge.from(2).to(4), 0),
						weighted(Edge.from(3).to(4), 0),
						weighted(Edge.from(4).to(2), 0)));
		final Weighted<Arborescence<Integer>> tree = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		// The reported score is not checked: 1e16 cancels the small weights while it is accumulated.
		assertEquals(Integer.valueOf(0), tree.val.parents.get(1));
		assertEquals(5, tree.val.parents.size());
	}

	@Test
	public void testPathsJoinPreviouslyProcessedBranches() {
		final WeightedGraph<Integer> graph = SparseWeightedGraph.from(
				ImmutableList.of(5, 4, 3, 2, 1, 0), ImmutableList.of(
						weighted(Edge.from(0).to(1), 4),
						weighted(Edge.from(1).to(2), 3),
						weighted(Edge.from(1).to(3), 2),
						weighted(Edge.from(3).to(4), 1),
						weighted(Edge.from(2).to(5), -2),
						weighted(Edge.from(4).to(4), 100),
						weighted(Edge.from(5).to(0), 100)));
		final Weighted<Arborescence<Integer>> tree = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		assertRootedTree(graph, 0, tree);
		assertEquals(8.0, tree.weight, DELTA);
	}

	@Test
	public void testLongPathFollowedByManySingletonPaths() {
		final int length = 1024;
		final List<Integer> nodes = new ArrayList<>();
		final List<Weighted<Edge<Integer>>> edges = new ArrayList<>();
		final ImmutableMap.Builder<Integer, Integer> parents = ImmutableMap.builder();
		for (int node = 0; node <= 2 * length; node++) nodes.add(node);
		for (int destination = 0; destination < nodes.size(); destination++) {
			if (destination == length) continue;
			final Integer source = nodes.get(destination < length ? destination + 1 : length);
			edges.add(weighted(Edge.from(source).to(nodes.get(destination)), 1.0));
			parents.put(nodes.get(destination), source);
		}
		final WeightedGraph<Integer> graph = SparseWeightedGraph.from(nodes, edges);
		final Weighted<Arborescence<Integer>> tree =
				ChuLiuEdmonds.getMaxArborescence(graph, nodes.get(length));
		assertEquals(parents.build(), tree.val.parents);
		assertEquals(2.0 * length, tree.weight, DELTA);
	}

	@Test
	public void testDisconnectedForestAndIsolatedNode() {
		final WeightedGraph<Integer> graph = SparseWeightedGraph.from(
				ImmutableList.of(4, 3, 2, 1, 0, 5), ImmutableList.of(
						weighted(Edge.from(0).to(1), 3),
						weighted(Edge.from(1).to(2), -1),
						weighted(Edge.from(3).to(4), 2)));
		final Weighted<Arborescence<Integer>> forest = ChuLiuEdmonds.getMaxArborescence(graph);
		assertEquals(ImmutableMap.of(1, 0, 2, 1, 4, 3), forest.val.parents);
		assertEquals(4.0, forest.weight, DELTA);
		assertEdgesSumToScore(graph, forest);
	}

	@Test
	public void testEmptyAndSingletonGraphs() {
		final WeightedGraph<Integer> empty = SparseWeightedGraph.from(
				ImmutableList.<Integer>of(), ImmutableList.<Weighted<Edge<Integer>>>of());
		assertEquals(weighted(Arborescence.<Integer>empty(), 0.0), ChuLiuEdmonds.getMaxArborescence(empty));
		final WeightedGraph<Integer> singleton = SparseWeightedGraph.from(
				ImmutableList.of(0), ImmutableList.of(weighted(Edge.from(0).to(0), 100)));
		assertRootedTree(singleton, 0, ChuLiuEdmonds.getMaxArborescence(singleton, 0));
	}

	@Test
	public void testDenseAndSparseGraphsAgainstExhaustiveOracle() {
		final Random random = new Random(20261101);
		for (int sample = 0; sample < 60; sample++) {
			final int size = 2 + sample % 5;
			final int root = random.nextInt(size);
			final double[][] weights = new double[size][size];
			final List<Integer> nodes = new ArrayList<>();
			final List<Weighted<Edge<Integer>>> edges = new ArrayList<>();
			for (int source = 0; source < size; source++) {
				nodes.add(source);
				for (int destination = 0; destination < size; destination++) {
					final double weight = source == root || sample % 2 == 0 || random.nextBoolean()
							? random.nextInt(11) - 5 : NINF;
					weights[source][destination] = weight;
					if (weight != NINF) edges.add(weighted(Edge.from(source).to(destination), weight));
				}
			}
			Collections.shuffle(nodes, random);
			final WeightedGraph<Integer> dense = DenseWeightedGraph.from(weights);
			final WeightedGraph<Integer> sparse = SparseWeightedGraph.from(nodes, edges);
			final double expected = enumerateArborescences(dense, root).stream()
					.mapToDouble(tree -> tree.weight).max().getAsDouble();
			for (WeightedGraph<Integer> graph : ImmutableList.of(dense, sparse)) {
				final Weighted<Arborescence<Integer>> actual = ChuLiuEdmonds.getMaxArborescence(graph, root);
				assertEquals("Sample " + sample + ", root " + root, expected, actual.weight, DELTA);
				assertRootedTree(graph, root, actual);
			}
		}
	}

	@Test
	public void testConstraintsAgainstExhaustiveOracle() {
		final Random random = new Random(481516);
		final Edge<Integer> required = Edge.from(0).to(1);
		final Edge<Integer> banned = Edge.from(2).to(3);
		for (int sample = 0; sample < 30; sample++) {
			final List<Weighted<Edge<Integer>>> edges = new ArrayList<>();
			for (int source = 0; source < 4; source++) {
				for (int destination = 1; destination < 4; destination++) {
					if (source != destination) {
						edges.add(weighted(Edge.from(source).to(destination), random.nextInt(11) - 5));
					}
				}
			}
			final WeightedGraph<Integer> graph = SparseWeightedGraph.from(edges);
			final double expected = enumerateArborescences(graph, 0).stream()
					.filter(tree -> tree.val.contains(required) && !tree.val.contains(banned))
					.mapToDouble(tree -> tree.weight).max().getAsDouble();
			final Weighted<Arborescence<Integer>> actual = ChuLiuEdmonds.getMaxArborescence(
					graph, ImmutableSet.of(required), ImmutableSet.of(banned));
			assertEquals("Constrained sample " + sample, expected, actual.weight, DELTA);
			assertTrue(actual.val.contains(required));
			assertFalse(actual.val.contains(banned));
			assertRootedTree(graph, 0, actual);
		}
	}

	@Test
	public void testGetMaxSpanningTree() {
		/*
		root    10
		(0) -------> (1) \
		 |  \       /  ^  \
		 |   \30   |   |20 \
		 |10  \    |10 |    \10
		 |     \   |  /      \
		 V  15  V  V /   20   V
		(3)<----- (2) -----> (4)
		  \-------^
		     40
		 */
		double[][] weights = {
				{NINF, 10, 30, 10, NINF},
				{NINF, NINF, 10, NINF, 10 },
				{NINF,  20, NINF,  7, 20 },
				{NINF, NINF, 40, NINF, NINF},
				{NINF, NINF, NINF, NINF, NINF},
		};
		final DenseWeightedGraph<Integer> graph = DenseWeightedGraph.from(weights);
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		/*
		root
		(0)           (1)
		 |             ^
		 |             |
		 |             |
		 |            /
		 V           /
		(3)       (2) ------> (4)
		  \-------^
		 */
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val.parents;
		assertEquals(2, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(2, maxBranching.get(4).intValue());
		assertEquals(90.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);
	}

	@Test
	public void testRequiredAndBannedEdges() {
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxArborescence(
				graph,
				ImmutableSet.of(Edge.from(0).to(1)),
				ImmutableSet.of(Edge.from(2).to(3)));
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val.parents;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(1, maxBranching.get(3).intValue());
		assertEquals(20.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);

	}

	@Test
	public void testRequiredAndBannedEdges2() {
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxArborescence(
				graph,
				ImmutableSet.of(Edge.from(0).to(3), Edge.from(3).to(1)),
				ImmutableSet.of(Edge.from(1).to(2))
		);
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val.parents;
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(18.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);

	}

	@Test
	public void testElevenNodeGraph() {
		// make a graph with a bunch of nested cycles so we can exercise the recursive part of the algorithm.
		final WeightedGraph<Integer> graph = SparseWeightedGraph.from(ImmutableList.of(
				weighted(Edge.from(0).to(8), 0),
				weighted(Edge.from(1).to(2), 10),
				weighted(Edge.from(1).to(4), 5),
				weighted(Edge.from(2).to(3), 9),
				weighted(Edge.from(3).to(1), 8),
				weighted(Edge.from(4).to(5), 9),
				weighted(Edge.from(5).to(6), 10),
				weighted(Edge.from(6).to(4), 8),
				weighted(Edge.from(6).to(7), 5),
				weighted(Edge.from(7).to(8), 10),
				weighted(Edge.from(8).to(2), 5),
				weighted(Edge.from(8).to(9), 8),
				weighted(Edge.from(8).to(10), 1),
				weighted(Edge.from(9).to(7), 9),
				weighted(Edge.from(10).to(3), 3)
		));
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxArborescence(graph, 0);

		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val.parents;
		assertEdgesSumToScore(graph, weightedSpanningTree);
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(8, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(1, maxBranching.get(4).intValue());
		assertEquals(4, maxBranching.get(5).intValue());
		assertEquals(5, maxBranching.get(6).intValue());
		assertEquals(9, maxBranching.get(7).intValue());
		assertEquals(0, maxBranching.get(8).intValue());
		assertEquals(8, maxBranching.get(9).intValue());
		assertEquals(8, maxBranching.get(10).intValue());
	}
}
