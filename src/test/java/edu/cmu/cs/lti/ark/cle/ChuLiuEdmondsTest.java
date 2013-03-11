package edu.cmu.cs.lti.ark.cle;

import org.junit.Test;

import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * @author sthomson@cs.cmu.edu
 */
public class ChuLiuEdmondsTest {
	@Test
	public void testGetMaxBranching() {
		double ninf = Double.NEGATIVE_INFINITY;
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
				{ ninf, 10, 30, 10, ninf },
				{ ninf, ninf, 10, ninf, 10 },
				{ ninf,  20, ninf,  7, 20 },
				{ ninf, ninf, 40, ninf, ninf },
				{ ninf, ninf, ninf, ninf, ninf},
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
	public void testGetMaxBranching2() {
		double ninf = Double.NEGATIVE_INFINITY;
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
				{ ninf, 10, ninf, 10, ninf },
				{ ninf, ninf, ninf, ninf, 10 },
				{ ninf,  ninf, ninf,  35, ninf },
				{ ninf, 20, 40, ninf, ninf },
				{ ninf, 20, 50, ninf, ninf},
		};
		final Weighted<Map<Integer, Integer>> weightedSpanningTree = ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		Map<Integer, Integer> maxBranching = weightedSpanningTree.val;
				// Print maximum branching per node.
		System.out.println("Maximum branching:");
		for (int to = 1; to <= maxBranching.size(); ++to)
			System.out.println(maxBranching.get(to) + " -> " + to);
		System.out.println(weightedSpanningTree.weight);

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
		assertEquals(0, maxBranching.get(1).intValue());
		assertEquals(4, maxBranching.get(2).intValue());
		assertEquals(2, maxBranching.get(3).intValue());
		assertEquals(1, maxBranching.get(4).intValue());
		assertEquals(105.0, weightedSpanningTree.weight);
	}
}
