package edu.cmu.cs.ark.cle;

import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.*;
import org.junit.Test;

import java.util.List;
import java.util.Map;
import java.util.Set;

import static edu.cmu.cs.ark.cle.ChuLiuEdmonds.SubsetOfSolutions;
import static edu.cmu.cs.ark.cle.EdgeQueueMap.EdgeQueue;
import static edu.cmu.cs.ark.cle.Weighted.weighted;
import static org.junit.Assert.*;

/**
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmondsTest {
	private final static double DELTA = 0.001;
	private final static double NINF = Double.NEGATIVE_INFINITY;
	private final static ImmutableSet<Edge<Integer>> empty = ImmutableSet.of();
	private final static WeightedGraph<Integer> graph = SparseWeightedGraph.from(ImmutableList.of(
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

	private <V> void assertEdgesSumToScore(WeightedGraph<V> originalEdgeWeights, Weighted<Arborescence<V>> bestTree) {
		final Map<V, V> parentsMap = bestTree.val.parents;
		double sumOfWeights = 0.0;
		for (V dest : parentsMap.keySet()) {
			final V source = parentsMap.get(dest);
			sumOfWeights += originalEdgeWeights.getWeightOf(source, dest);
		}
		assertEquals(sumOfWeights, bestTree.weight, DELTA);
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
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(graph, 0);
		printTree(weightedSpanningTree);
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
	public void testGetKBest() {
		final List<Weighted<Arborescence<Integer>>> weightedSpanningTrees =
				ChuLiuEdmonds.getKBestSpanningTrees(graph, 0, 4);

		Weighted<Arborescence<Integer>> weightedSpanningTree = weightedSpanningTrees.get(0);
		Map<Integer, Integer> maxBranching = weightedSpanningTree.val.parents;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(21.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(1);
		maxBranching = weightedSpanningTree.val.parents;
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(21.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(2);
		maxBranching = weightedSpanningTree.val.parents;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(1, maxBranching.get(3).intValue());
		assertEquals(20.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(3);
		maxBranching = weightedSpanningTree.val.parents;
		assertEquals(2, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(19.0, weightedSpanningTrees.get(3).weight, DELTA);
		assertEdgesSumToScore(graph, weightedSpanningTree);
	}

	@Test
	public void testGetLotsOfKBest() {
		final int k = 100;
		final List<Weighted<Arborescence<Integer>>> kBestSpanningTrees = ChuLiuEdmonds.getKBestSpanningTrees(graph, 0, k);
		final int size = kBestSpanningTrees.size();
		// make sure there are no more than k of them
		assertTrue(size <= k);
		// make sure they are in descending order
		for (int i = 0; i + 1 < size; i++) {
			assertTrue(kBestSpanningTrees.get(i).weight >= kBestSpanningTrees.get(i+1).weight);
		}
		// make sure they're all unique
		final Set<String> kBestStrings = Sets.newHashSet();
		for (Weighted<Arborescence<Integer>> spanningTree : kBestSpanningTrees) {
			kBestStrings.add(showTree(spanningTree));
		}
		assertEquals(size, kBestStrings.size());
	}

	@Test
	public void testSeekDoesntReturnAncestor() {
		final Weighted<Arborescence<Integer>> bestArborescence = ChuLiuEdmonds.getMaxSpanningTree(graph, 0);
		final ExclusiveEdge<Integer> maxInEdge = ExclusiveEdge.of(Edge.from(1).to(2), 11.0);
		final EdgeQueue<Integer> edgeQueue = EdgeQueue.create(maxInEdge.edge.destination, Partition.singletons(graph.getNodes()));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(0).to(2), 1.0));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(3).to(2), 8.0));
		final Optional<ExclusiveEdge<Integer>> nextBestEdge = ChuLiuEdmonds.seek(maxInEdge, bestArborescence.val, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		// 3 -> 2 is an ancestor in bestArborescence, so seek should not return it
		assertNotEquals(Edge.from(3).to(2), nextBestEdge.get().edge);
		assertEquals(Edge.from(0).to(2), nextBestEdge.get().edge);
	}

	@Test
	public void testSeek() {
		final Arborescence<Integer> best = Arborescence.of(ImmutableMap.of(
				2, 0,
				1, 2,
				3, 2
		));
		final ExclusiveEdge<Integer> maxInEdge = ExclusiveEdge.of(Edge.from(2).to(1), 10.0);
		final EdgeQueue<Integer> edgeQueue = EdgeQueue.create(maxInEdge.edge.destination, Partition.singletons(graph.getNodes()));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), 5.0));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), 9.0));
		final Optional<ExclusiveEdge<Integer>> nextBestEdge = ChuLiuEdmonds.seek(maxInEdge, best, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		assertEquals(Edge.from(3).to(1), nextBestEdge.get().edge);
		assertEquals(9.0, nextBestEdge.get().weight, DELTA);
	}

	@Test
	public void testNext() {
		// get the best tree A(1)
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxSpanningTree(graph, 0);
		final Optional<SubsetOfSolutions<Integer>> oItem =
				ChuLiuEdmonds.getNextBest(graph, 0, empty, empty, best);
		assertTrue(oItem.isPresent());
		final SubsetOfSolutions<Integer> item = oItem.get();
		assertEquals(Edge.from(0).to(1), item.edgeToBan);
		assertEquals(0.0, item.bestArborescence.weight - item.weightOfNextBest, DELTA);
	}

	@Test
	public void testNextWithRequiredEdges() {
		// get the best tree A(1)
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxSpanningTree(graph, 0);
		final Optional<SubsetOfSolutions<Integer>> oItem =
				ChuLiuEdmonds.getNextBest(graph, 0, ImmutableSet.of(Edge.from(0).to(1)), empty, best);
		assertTrue(oItem.isPresent());
		final SubsetOfSolutions<Integer> item = oItem.get();
		assertEquals(Edge.from(2).to(3), item.edgeToBan);
		assertEquals(1.0, item.bestArborescence.weight - item.weightOfNextBest, DELTA);
	}

	@Test
	public void testNextReturnsAbsentWhenTreesAreExhausted() {
		// get the best tree A(1)
		final WeightedGraph<Integer> graph = DenseWeightedGraph.from(new double[][]{
				{NINF, 1.0},
				{NINF, NINF}
		});
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxSpanningTree(graph, 0);
		Optional<SubsetOfSolutions<Integer>> pair = ChuLiuEdmonds.getNextBest(graph, 0, empty, empty, best);
		assertFalse(pair.isPresent());
	}

	@Test
	public void testRequiredAndBannedEdges() {
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(
				graph,
				0,
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
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(
				graph,
				0,
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

	private <V> void printTree(Weighted<Arborescence<V>> weightedSpanningTree) {
		System.out.println(showTree(weightedSpanningTree));
	}

	private <V> String showTree(Weighted<Arborescence<V>> weightedSpanningTree) {
		List<String> lines = Lists.newArrayList();
		final Map<V, V> maxBranching = weightedSpanningTree.val.parents;
		for (V to : maxBranching.keySet()) {
			lines.add(maxBranching.get(to) + " -> " + to);
		}
		lines.add(Double.toString(weightedSpanningTree.weight));
		return Joiner.on("\n").join(lines);
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
		final Weighted<Arborescence<Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(graph, 0);
		printTree(weightedSpanningTree);

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
