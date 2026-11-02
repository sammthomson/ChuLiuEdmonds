package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import edu.cmu.cs.ark.cle.ds.Partition;
import edu.cmu.cs.ark.cle.graph.Edge;
import edu.cmu.cs.ark.cle.graph.SparseWeightedGraph;
import edu.cmu.cs.ark.cle.graph.WeightedGraph;
import edu.cmu.cs.ark.cle.util.Weighted;
import org.junit.Test;

import java.util.*;

import static edu.cmu.cs.ark.cle.ChuLiuEdmondsTest.*;
import static edu.cmu.cs.ark.cle.util.Weighted.weighted;
import static org.junit.Assert.*;

public class KBestArborescencesTest {
	private final static ImmutableSet<Edge<Integer>> empty = ImmutableSet.of();
	// tied for first, can appear in either order
	private final static Weighted<Arborescence<Integer>> bestA = weighted(Arborescence.of(ImmutableMap.of(
			1, 0,
			2, 1,
			3, 2
	)), 21);
	private final static Weighted<Arborescence<Integer>> bestB = weighted(Arborescence.of(ImmutableMap.of(
			1, 3,
			2, 1,
			3, 0
	)), 21);
	final static Set<Weighted<Arborescence<Integer>>> expectedFirstAndSecond = ImmutableSet.of(bestA, bestB);

	private static void assertRankingMatchesOracle(WeightedGraph<Integer> graph, int root) {
		final List<Weighted<Arborescence<Integer>>> expected = enumerateArborescences(graph, root);
		expected.sort(Comparator.comparingDouble((Weighted<Arborescence<Integer>> tree) -> tree.weight).reversed());
		final List<Weighted<Arborescence<Integer>>> actual =
				KBestArborescences.getKBestArborescences(graph, root, expected.size() + 1);
		assertEquals(expected.size(), actual.size());
		assertEquals(new HashSet<>(expected), new HashSet<>(actual));
		for (int index = 0; index < expected.size(); index++) {
			assertEquals("Rank " + index, expected.get(index).weight, actual.get(index).weight, DELTA);
			assertRootedTree(graph, root, actual.get(index));
		}
	}

	@Test
	public void testNestedCyclesRankingAgainstExhaustiveOracle() {
		assertRankingMatchesOracle(nestedCyclesGraph(), 0);
	}

	@Test
	public void testTiedRankingsAgainstExhaustiveOracle() {
		final Random random = new Random(271828);
		for (int sample = 0; sample < 20; sample++) {
			final List<Weighted<Edge<Integer>>> edges = new ArrayList<>();
			for (int source = 0; source < 4; source++) {
				for (int destination = 0; destination < 4; destination++) {
					if (source != destination) {
						edges.add(weighted(Edge.from(source).to(destination), random.nextInt(5) - 2));
					}
				}
			}
			assertRankingMatchesOracle(SparseWeightedGraph.from(edges), sample % 4);
		}
	}

	@Test
	public void testGetKBestArborescences() {
		final List<Weighted<Arborescence<Integer>>> weightedSpanningTrees =
				KBestArborescences.getKBestArborescences(graph, 0, 4);

		assertEquals(4, ImmutableSet.copyOf(weightedSpanningTrees).size());

		Weighted<Arborescence<Integer>> weightedSpanningTree = weightedSpanningTrees.get(0);
		assertTrue(expectedFirstAndSecond.contains(weightedSpanningTree));
		assertEdgesSumToScore(graph, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(1);
		assertTrue(expectedFirstAndSecond.contains(weightedSpanningTree));
		assertEdgesSumToScore(graph, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(2);
		final Arborescence<Integer> expectedThird = Arborescence.of(ImmutableMap.of(
				1, 0,
				2, 1,
				3, 1));
		assertEquals(weighted(expectedThird, 20), weightedSpanningTree);
		assertEdgesSumToScore(graph, weightedSpanningTree);

		weightedSpanningTree = weightedSpanningTrees.get(3);
		final Arborescence<Integer> expectedFourth = Arborescence.of(ImmutableMap.of(
				1, 2,
				2, 3,
				3, 0));
		assertEquals(weighted(expectedFourth, 19.0), weightedSpanningTree);
		assertEdgesSumToScore(graph, weightedSpanningTree);
	}

	@Test
	public void testGetLotsOfKBest() {
		final int k = 100;
		final List<Weighted<Arborescence<Integer>>> kBestSpanningTrees =
				KBestArborescences.getKBestArborescences(graph, 0, k);
		final int size = kBestSpanningTrees.size();
		// make sure there are no more than k of them
		assertTrue(size <= k);
		// make sure they are in descending order
		for (int i = 0; i + 1 < size; i++) {
			assertTrue(kBestSpanningTrees.get(i).weight >= kBestSpanningTrees.get(i+1).weight);
		}
		for (Weighted<Arborescence<Integer>> spanningTree : kBestSpanningTrees) {
			assertEdgesSumToScore(graph, spanningTree);
		}
		// make sure they're all unique
		assertEquals(size, ImmutableSet.copyOf(kBestSpanningTrees).size());
	}

	@Test
	public void testSeekDoesntReturnAncestor() {
		final Weighted<Arborescence<Integer>> bestArborescence = bestA;
		final ExclusiveEdge<Integer> maxInEdge = ExclusiveEdge.of(Edge.from(1).to(2), 11.0);
		final EdgeQueueMap.EdgeQueue<Integer> edgeQueue =
				EdgeQueueMap.EdgeQueue.create(maxInEdge.edge.destination, Partition.singletons(graph.getNodes()));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(0).to(2), 1.0));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(3).to(2), 8.0));
		final Optional<ExclusiveEdge<Integer>> nextBestEdge =
				KBestArborescences.seek(maxInEdge, bestArborescence.val, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		// 3 -> 2 is an ancestor in bestArborescence, so seek should not return it
		assertNotEquals(Edge.from(3).to(2), nextBestEdge.get().edge);
		assertEquals(Edge.from(0).to(2), nextBestEdge.get().edge);
	}

	@Test
	public void testSeekPreservesAncestorEdgesForCycleContraction() {
		for (boolean hasAlternative : new boolean[] {false, true}) {
			final ExclusiveEdge<Integer> maxInEdge = ExclusiveEdge.of(Edge.from(1).to(2), 11.0);
			final EdgeQueueMap.EdgeQueue<Integer> queue =
					EdgeQueueMap.EdgeQueue.create(2, Partition.singletons(graph.getNodes()));
			final ExclusiveEdge<Integer> ancestorEdge = ExclusiveEdge.of(Edge.from(3).to(2), 8.0);
			queue.addEdge(ancestorEdge);
			if (hasAlternative) queue.addEdge(ExclusiveEdge.of(Edge.from(0).to(2), 1.0));
			assertEquals(hasAlternative, KBestArborescences.seek(maxInEdge, bestA.val, queue).isPresent());
			assertEquals(ancestorEdge.edge, queue.popBestEdge().get().edge);
		}
	}

	@Test
	public void testSeek() {
		final Arborescence<Integer> best = Arborescence.of(ImmutableMap.of(
				2, 0,
				1, 2,
				3, 2
		));
		final ExclusiveEdge<Integer> maxInEdge = ExclusiveEdge.of(Edge.from(2).to(1), 10.0);
		final EdgeQueueMap.EdgeQueue<Integer> edgeQueue =
				EdgeQueueMap.EdgeQueue.create(maxInEdge.edge.destination, Partition.singletons(graph.getNodes()));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), 5.0));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), 9.0));
		final Optional<ExclusiveEdge<Integer>> nextBestEdge = KBestArborescences.seek(maxInEdge, best, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		assertEquals(Edge.from(3).to(1), nextBestEdge.get().edge);
		assertEquals(9.0, nextBestEdge.get().weight, DELTA);
	}

	@Test
	public void testNext() {
		final Optional<Weighted<KBestArborescences.SubsetOfSolutions<Integer>>> oItem =
				KBestArborescences.scoreSubsetOfSolutions(graph, empty, empty, bestA);
		assertTrue(oItem.isPresent());
		final KBestArborescences.SubsetOfSolutions<Integer> item = oItem.get().val;
		assertEquals(Edge.from(0).to(1), item.edgeToBan);
		assertEquals(0.0, item.bestArborescence.weight - oItem.get().weight, DELTA);
	}

	@Test
	public void testNextWithRequiredEdges() {
		final Optional<Weighted<KBestArborescences.SubsetOfSolutions<Integer>>> oItem =
				KBestArborescences.scoreSubsetOfSolutions(graph, ImmutableSet.of(Edge.from(0).to(1)), empty, bestA);
		assertTrue(oItem.isPresent());
		final KBestArborescences.SubsetOfSolutions<Integer> item = oItem.get().val;
		assertEquals(Edge.from(2).to(3), item.edgeToBan);
		assertEquals(1.0, item.bestArborescence.weight - oItem.get().weight, DELTA);
	}

	@Test
	public void testNextReturnsAbsentWhenTreesAreExhausted() {
		final WeightedGraph<Integer> oneTreeGraph = SparseWeightedGraph.from(
				ImmutableSet.of(weighted(Edge.from(0).to(1), 1.0))
		);
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxArborescence(oneTreeGraph, 0);
		final Optional<Weighted<KBestArborescences.SubsetOfSolutions<Integer>>> pair =
				KBestArborescences.scoreSubsetOfSolutions(oneTreeGraph, empty, empty, best);
		assertFalse(pair.isPresent());
	}
}
