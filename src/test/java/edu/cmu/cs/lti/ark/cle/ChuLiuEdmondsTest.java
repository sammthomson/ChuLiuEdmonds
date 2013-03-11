package edu.cmu.cs.lti.ark.cle;

import org.junit.Test;

import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmondsTest {
	private double NINF = Double.NEGATIVE_INFINITY;

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
		Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
				// Print maximum branching per node.
		System.out.println("Maximum branching:");
		for (int to = 1; to <= maxBranching.size(); ++to)
			System.out.println(maxBranching.get(to) + " -> " + to);
		System.out.println(weightedSpanningTree.weight);
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
		assertEquals(2, maxBranching.get(1).intValue());
		assertEquals(3, maxBranching.get(2).intValue());
		assertEquals(0, maxBranching.get(3).intValue());
		assertEquals(2, maxBranching.get(4).intValue());
		assertEquals(90.0, weightedSpanningTree.weight);
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
		// Print maximum branching per node.
		System.out.println("Maximum branchings:");
		Map<Integer, Integer> maxBranching = weightedSpanningTrees.get(0).val;
		for (Weighted<Map<Integer, Integer>> weightedSpanningTree : weightedSpanningTrees) {
			maxBranching = weightedSpanningTree.val;
			for (int to = 1; to <= maxBranching.size(); ++to) {
				System.out.println(maxBranching.get(to) + " -> " + to);
			}
			System.out.println(weightedSpanningTree.weight);
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
		maxBranching = weightedSpanningTrees.get(0).val;
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(4, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(1, maxBranching.get(4).intValue());
		assertEquals(105.0, weightedSpanningTrees.get(0).weight);
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
		assertEquals(30.0, weightedSpanningTrees.get(1).weight);

	}
}
