package edu.cmu.cs.ark.cle;

import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.Test;

import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.*;

/**
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmondsTest {
	private final static double DELTA = 0.001;
	private final static double NINF = Double.NEGATIVE_INFINITY;
	private final double[][] weights = new double[][] {
			{NINF, 5, 1, 1},
			{NINF, NINF, 11, 4},
			{NINF,  10, NINF,  5},
			{NINF, 9, 8, NINF},
	};

	private void assertEdgesSumToScore(double[][] originalEdgeWeights, Weighted<Map<Integer,Integer>> parentsMap) {
		double sumOfWeights = 0.0;
		for (int to = 1; to <= parentsMap.val.size(); ++to) {
			final int from = parentsMap.val.get(to);
			sumOfWeights += originalEdgeWeights[from][to];
		}
		assertEquals(sumOfWeights, parentsMap.weight, DELTA);
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
		final Weighted<Map<Integer, Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
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
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEquals(2, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(2, maxBranching.get(4).intValue());
		assertEquals(90.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);
	}

	@Test
	public void testGetKBest() {
		final List<Weighted<Map<Integer, Integer>>> weightedSpanningTrees =
				ChuLiuEdmonds.getKBestSpanningTrees(weights, 0, 4);

		Weighted<Map<Integer, Integer>> weightedSpanningTree = weightedSpanningTrees.get(0);
		Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(21.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(1);
		maxBranching = weightedSpanningTree.val;
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(21.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(2);
		maxBranching = weightedSpanningTree.val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(1, maxBranching.get(3).intValue());
		assertEquals(20.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(3);
		maxBranching = weightedSpanningTree.val;
		assertEquals(2, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(19.0, weightedSpanningTrees.get(3).weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);
	}

	@Test
	public void testGetLotsOfKBest() {
		final int k = 100;
		final List<Weighted<Map<Integer, Integer>>> kBestSpanningTrees = ChuLiuEdmonds.getKBestSpanningTrees(weights, 0, k);
		final int size = kBestSpanningTrees.size();
		// make sure there are no more than k of them
		assertTrue(size <= k);
		// make sure they are in descending order
		for (int i = 0; i + 1 < size; i++) {
			assertTrue(kBestSpanningTrees.get(i).weight >= kBestSpanningTrees.get(i+1).weight);
		}
		// make sure they're all unique
		final Set<String> kBestStrings = Sets.newHashSet();
		for (Weighted<Map<Integer, Integer>> spanningTree : kBestSpanningTrees) {
			kBestStrings.add(showTree(spanningTree));
		}
		assertEquals(size, kBestStrings.size());
	}

	@Test
	public void testSeekDoesntReturnAncestor() {
		final Weighted<Map<Integer, Integer>> bestArborescence = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		final ExclusiveEdge maxInEdge = new ExclusiveEdge(new Edge(1, 2), ImmutableList.<Edge>of(), 11.0);
		final EdgeQueueMap.EdgeQueue edgeQueue = new EdgeQueueMap.EdgeQueue(maxInEdge.edge.destination, new Partition(4));
		edgeQueue.addEdge(new ExclusiveEdge(new Edge(0, 2), ImmutableList.<Edge>of(), 1.0));
		edgeQueue.addEdge(new ExclusiveEdge(new Edge(3, 2), ImmutableList.<Edge>of(), 8.0));
		final Optional<ExclusiveEdge> nextBestEdge = ChuLiuEdmonds.seek(maxInEdge, bestArborescence.val, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		// 3 -> 2 is an ancestor in bestArborescence, so seek should not return it
		assertNotEquals(new Edge(3, 2), nextBestEdge.get().edge);
		assertEquals(new Edge(0, 2), nextBestEdge.get().edge);
	}

	@Test
	public void testSeek() {
		final Map<Integer, Integer> best = ImmutableMap.of(
				2, 0,
				1, 2,
				3, 2
		);
		final ExclusiveEdge maxInEdge = new ExclusiveEdge(new Edge(2, 1), ImmutableList.<Edge>of(), 10.0);
		final EdgeQueueMap.EdgeQueue edgeQueue = new EdgeQueueMap.EdgeQueue(maxInEdge.edge.destination, new Partition(4));
		edgeQueue.addEdge(new ExclusiveEdge(new Edge(0, 1), ImmutableList.<Edge>of(), 5.0));
		edgeQueue.addEdge(new ExclusiveEdge(new Edge(3, 1), ImmutableList.<Edge>of(), 9.0));
		final Optional<ExclusiveEdge> nextBestEdge = ChuLiuEdmonds.seek(maxInEdge, best, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		assertEquals(new Edge(3, 1), nextBestEdge.get().edge);
		assertEquals(9.0, nextBestEdge.get().weight, DELTA);
	}

	@Test
	public void testNext() {
		// get the best tree A(1)
		final Weighted<Map<Integer, Integer>> best = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		Optional<Pair<Edge, Double>> oPair =
				ChuLiuEdmonds.next(weights, 0, ImmutableList.<Edge>of(), ImmutableList.<Edge>of(), best);
		assertTrue(oPair.isPresent());
		final Pair<Edge, Double> pair = oPair.get();
		assertEquals(new Edge(0, 1), pair.first);
		assertEquals(0.0, pair.second, DELTA);
	}

	@Test
	public void testNextWithRequiredEdges() {
		// get the best tree A(1)
		final Weighted<Map<Integer, Integer>> best = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		Optional<Pair<Edge, Double>> oPair =
				ChuLiuEdmonds.next(weights, 0, ImmutableList.of(new Edge(0, 1)), ImmutableList.<Edge>of(), best);
		assertTrue(oPair.isPresent());
		final Pair<Edge, Double> pair = oPair.get();
		assertEquals(new Edge(2, 3), pair.first);
		assertEquals(1.0, pair.second, DELTA);
	}

	@Test
	public void testNextReturnsAbsentWhenTreesAreExhausted() {
		final double[][] weights = {
				{NINF, 1.0},
				{NINF, NINF}
		};
		// get the best tree A(1)
		final Weighted<Map<Integer, Integer>> best = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		Optional<Pair<Edge, Double>> pair =
				ChuLiuEdmonds.next(weights, 0, ImmutableList.<Edge>of(), ImmutableList.<Edge>of(), best);
		assertFalse(pair.isPresent());
	}

	@Test
	public void testRequiredAndBannedEdges() {
		final Weighted<Map<Integer, Integer>> weightedSpanningTree =
				ChuLiuEdmonds.getMaxSpanningTree(weights, 0, ImmutableList.of(new Edge(0, 1)), ImmutableList.of(new Edge(2, 3)));
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(1, maxBranching.get(3).intValue());
		assertEquals(20.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

	}

	@Test
	public void testRequiredAndBannedEdges2() {
		final Weighted<Map<Integer, Integer>> weightedSpanningTree =
				ChuLiuEdmonds.getMaxSpanningTree(weights, 0, ImmutableList.of(new Edge(0, 3), new Edge(3, 1)), ImmutableList.of(new Edge(1, 2)));
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(18.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

	}

	@Test
	public void testGetDiverseKBestMaxSpanningTrees() {
		/*
		root    10
		(0) -------> (1) <
		 |     ---^   |  \
		 |10  /        |10 \
		 |   /20        \   \ 20
		 |  /            \   \
		 V /  35           V  |
		(3)<----- (2) <----- (4)
		  \-------^     50
		     40
		 */
		double[][] weights = {
				{NINF, 10, NINF, 10, NINF},
				{NINF, NINF, NINF, NINF, 10 },
				{NINF, NINF, NINF,  35, NINF},
				{NINF, 20, 40, NINF, NINF},
				{NINF, 20, 50, NINF, NINF},
		};
		final List<Weighted<Map<Integer, Integer>>> weightedSpanningTrees =
				ChuLiuEdmonds.getDiverseKBestSpanningTrees(weights, 0, 2, 50.0);
		for (Weighted<Map<Integer, Integer>> weightedSpanningTree : weightedSpanningTrees) {
			printTree(weightedSpanningTree);
		}
		/*
		root    10
		(0) -------> (1)
		               |
		               |10
		                \
		                 \
		      35           V
		(3)<----- (2) <----- (4)
		                50
		 */
		final Weighted<Map<Integer, Integer>> weightedSpanningTree = weightedSpanningTrees.get(0);
		Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(4, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(1, maxBranching.get(4).intValue());
		assertEquals(105.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);
		/*
		root
		(0)         (1)
		 |     ---^   |
		 |10  /        |-40 (10 - 50)
		 |   /20        \
		 |  /            \
		 V /               V
		(3)       (2)       (4)
		  \-------^
		     40
		 */
		maxBranching = weightedSpanningTrees.get(1).val;
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(1, maxBranching.get(4).intValue());
		assertEquals(30.0, weightedSpanningTrees.get(1).weight, DELTA);

	}

	private void printTree(Weighted<Map<Integer, Integer>> weightedSpanningTree) {
		System.out.println(showTree(weightedSpanningTree));
	}

	private String showTree(Weighted<Map<Integer, Integer>> weightedSpanningTree) {
		List<String> lines = Lists.newArrayList();
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		for (int to = 1; to <= maxBranching.size(); to++) {
			lines.add(maxBranching.get(to) + " -> " + to);
		}
		lines.add(Double.toString(weightedSpanningTree.weight));
		return Joiner.on("\n").join(lines);
	}

	@Test
	public void testElevenNodeGraph() {
		// make a graph with a bunch of nested cycles so we can exercise the recursive part of the algorithm.
		double[][] weights = {
				{NINF, NINF, NINF, NINF, NINF, NINF, NINF, NINF, 0, NINF, NINF},
				{NINF, NINF, 10,   NINF, 5,    NINF, NINF, NINF, NINF, NINF, NINF}, //  \        \          \
				{NINF, NINF, NINF, 9,    NINF, NINF, NINF, NINF, NINF, NINF, NINF}, //  } cycle   \          \
				{NINF, 8,    NINF, NINF, NINF, NINF, NINF, NINF, NINF, NINF, NINF}, // /           \          \
				{NINF, NINF, NINF, NINF, NINF, 9,    NINF, NINF, NINF, NINF, NINF}, //  \           \          \
				{NINF, NINF, NINF, NINF, NINF, NINF, 10,   NINF, NINF, NINF, NINF}, //  } cycle     } cycle     \
				{NINF, NINF, NINF, NINF, 8,    NINF, NINF, 5,    NINF, NINF, NINF}, // /           /            } cycle
				{NINF, NINF, NINF, NINF, NINF, NINF, NINF, NINF, 10,   NINF, NINF}, //  \         /            /
				{NINF, NINF, 5   , NINF, NINF, NINF, NINF, NINF, NINF, 8,    1,  }, //  } cycle  /            /
				{NINF, NINF, NINF, NINF, NINF, NINF, NINF, 9,    NINF, NINF, NINF}, // /        /            /
				{NINF, NINF, NINF, 3   , NINF, NINF, NINF, NINF, NINF, NINF, NINF}, //                      /
		};
		final Weighted<Map<Integer, Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		printTree(weightedSpanningTree);

		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEdgesSumToScore(weights, weightedSpanningTree);
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
