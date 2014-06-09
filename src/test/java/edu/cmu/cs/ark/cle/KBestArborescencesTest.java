package edu.cmu.cs.ark.cle;

import com.google.common.base.Optional;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import org.junit.Test;

import java.util.List;
import java.util.Map;
import java.util.Set;

import static edu.cmu.cs.ark.cle.ChuLiuEdmondsTest.*;
import static org.junit.Assert.*;

public class KBestArborescencesTest {
	private final static ImmutableSet<Edge<Integer>> empty = ImmutableSet.of();

	@Test
	public void testGetKBestArborescences() {
		final List<Weighted<Arborescence<Integer>>> weightedSpanningTrees =
				KBestArborescences.getKBestArborescences(graph, 0, 4);

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
		final List<Weighted<Arborescence<Integer>>> kBestSpanningTrees = KBestArborescences.getKBestArborescences(graph, 0, k);
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
		final Weighted<Arborescence<Integer>> bestArborescence = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		final ExclusiveEdge<Integer> maxInEdge = ExclusiveEdge.of(Edge.from(1).to(2), 11.0);
		final EdgeQueueMap.EdgeQueue<Integer> edgeQueue = EdgeQueueMap.EdgeQueue.create(maxInEdge.edge.destination, Partition.singletons(graph.getNodes()));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(0).to(2), 1.0));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(3).to(2), 8.0));
		final Optional<ExclusiveEdge<Integer>> nextBestEdge = KBestArborescences.seek(maxInEdge, bestArborescence.val, edgeQueue);
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
		final EdgeQueueMap.EdgeQueue<Integer> edgeQueue = EdgeQueueMap.EdgeQueue.create(maxInEdge.edge.destination, Partition.singletons(graph.getNodes()));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), 5.0));
		edgeQueue.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), 9.0));
		final Optional<ExclusiveEdge<Integer>> nextBestEdge = KBestArborescences.seek(maxInEdge, best, edgeQueue);
		assertTrue(nextBestEdge.isPresent());
		assertEquals(Edge.from(3).to(1), nextBestEdge.get().edge);
		assertEquals(9.0, nextBestEdge.get().weight, DELTA);
	}

	@Test
	public void testNext() {
		// get the best tree A(1)
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		final Optional<Weighted<KBestArborescences.SubsetOfSolutions<Integer>>> oItem =
				KBestArborescences.scoreSubsetOfSolutions(graph, empty, empty, best);
		assertTrue(oItem.isPresent());
		final KBestArborescences.SubsetOfSolutions<Integer> item = oItem.get().val;
		assertEquals(Edge.from(0).to(1), item.edgeToBan);
		assertEquals(0.0, item.bestArborescence.weight - oItem.get().weight, DELTA);
	}

	@Test
	public void testNextWithRequiredEdges() {
		// get the best tree A(1)
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		final Optional<Weighted<KBestArborescences.SubsetOfSolutions<Integer>>> oItem =
				KBestArborescences.scoreSubsetOfSolutions(graph, ImmutableSet.of(Edge.from(0).to(1)), empty, best);
		assertTrue(oItem.isPresent());
		final KBestArborescences.SubsetOfSolutions<Integer> item = oItem.get().val;
		assertEquals(Edge.from(2).to(3), item.edgeToBan);
		assertEquals(1.0, item.bestArborescence.weight - oItem.get().weight, DELTA);
	}

	@Test
	public void testNextReturnsAbsentWhenTreesAreExhausted() {
		// get the best tree A(1)
		final WeightedGraph<Integer> graph = DenseWeightedGraph.from(new double[][]{
				{NINF, 1.0},
				{NINF, NINF}
		});
		final Weighted<Arborescence<Integer>> best = ChuLiuEdmonds.getMaxArborescence(graph, 0);
		Optional<Weighted<KBestArborescences.SubsetOfSolutions<Integer>>> pair =
				KBestArborescences.scoreSubsetOfSolutions(graph, empty, empty, best);
		assertFalse(pair.isPresent());
	}
}