package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableList;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmondsTest {
	private final static double DELTA = 0.001;
	private final static double NINF = Double.NEGATIVE_INFINITY;

	private void assertEdgesSumToScore(double[][] originalEdgeWeights, Weighted<Map<Integer,Integer>> parentsMap) {
		double sumOfWeights = 0.0;
		for (int to = 1; to <= parentsMap.val.size(); ++to) {
			final int from = parentsMap.val.get(to);
			sumOfWeights += originalEdgeWeights[from][to];
		}
		assertEquals(sumOfWeights, parentsMap.weight, 0.0001);
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
		double[][] weights = {
				{NINF, 5, 1, 1},
				{NINF, NINF, 11, 4},
				{NINF,  10, NINF,  5},
				{NINF, 9, 8, NINF},
		};

		final List<Weighted<Map<Integer, Integer>>> weightedSpanningTrees =
				ChuLiuEdmonds.getKBestSpanningTrees(weights, 0, 4);

		final Weighted<Map<Integer, Integer>> weightedSpanningTree = weightedSpanningTrees.get(0);
		Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(21.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

		maxBranching = weightedSpanningTrees.get(1).val;
		assertEquals(3, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(21.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

		maxBranching = weightedSpanningTrees.get(2).val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(1, maxBranching.get(2).intValue());
		assertEquals(1, maxBranching.get(3).intValue());
		assertEquals(20.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);

		maxBranching = weightedSpanningTrees.get(3).val;
		assertEquals(2, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(19.0, weightedSpanningTree.weight, DELTA);
		assertEdgesSumToScore(weights, weightedSpanningTree);
	}

	@Test
	public void testRequiredAndBannedEdges() {
		double[][] weights = {
				{NINF, 5, 1, 1},
				{NINF, NINF, 11, 4},
				{NINF,  10, NINF,  5},
				{NINF, 9, 8, NINF},
		};

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
		double[][] weights = {
				{NINF, 5, 1, 1},
				{NINF, NINF, 11, 4},
				{NINF,  10, NINF,  5},
				{NINF, 9, 8, NINF},
		};

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
		final Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
		for (int to = 1; to <= maxBranching.size(); to++) {
			System.out.println(maxBranching.get(to) + " -> " + to);
		}
		System.out.println(weightedSpanningTree.weight);
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
