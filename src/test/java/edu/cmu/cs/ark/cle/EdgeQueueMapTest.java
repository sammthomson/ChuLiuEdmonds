package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import edu.cmu.cs.ark.cle.ds.Partition;
import edu.cmu.cs.ark.cle.graph.Edge;
import edu.cmu.cs.ark.cle.util.Pair;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static edu.cmu.cs.ark.cle.EdgeQueueMap.EdgeQueue;
import static edu.cmu.cs.ark.cle.util.Weighted.weighted;
import static org.junit.Assert.*;

public class EdgeQueueMapTest {
	@Test
	public void testMergeAdjustsWeightsAndHistoriesWithoutMutatingSelectedEdges() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2, 3));
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		final Edge<Integer> previous = Edge.from(3).to(1);
		final ExclusiveEdge<Integer> original = ExclusiveEdge.of(
				Edge.from(0).to(1), ImmutableList.of(previous), 20.0);
		final EdgeQueue<Integer> first = EdgeQueue.create(1, partition);
		final EdgeQueue<Integer> second = EdgeQueue.create(2, partition);
		first.addEdge(original);
		first.addEdge(ExclusiveEdge.of(Edge.from(2).to(1), 100.0));
		second.addEdge(ExclusiveEdge.of(Edge.from(3).to(2), 19.0));

		final Edge<Integer> replaceFirst = Edge.from(2).to(1);
		final Edge<Integer> replaceSecond = Edge.from(1).to(2);
		final int component = partition.merge(1, 2);
		final EdgeQueue<Integer> merged = queues.merge(component, ImmutableList.of(
				Pair.of(first, weighted(replaceFirst, 8.0)),
				Pair.of(second, weighted(replaceSecond, 4.0))));

		assertFalse(first.popBestEdge().isPresent());
		assertFalse(second.popBestEdge().isPresent());
		final ExclusiveEdge<Integer> best = merged.popBestEdge().get();
		assertEquals(Edge.from(3).to(2), best.edge);
		assertEquals(15.0, best.weight, 0.0);
		assertEquals(ImmutableList.of(replaceSecond), best.excluded);
		final ExclusiveEdge<Integer> next = merged.popBestEdge().get();
		assertEquals(original.edge, next.edge);
		assertEquals(12.0, next.weight, 0.0);
		assertEquals(ImmutableList.of(previous, replaceFirst), next.excluded);
		assertEquals(ImmutableList.of(previous), original.excluded);
		assertFalse(merged.popBestEdge().isPresent());

		merged.addEdge(next);
		assertEquals(next.weight, merged.popBestEdge().get().weight, 0.0);
		merged.addEdge(best);
		assertEquals(best.excluded, merged.popBestEdge().get().excluded);
	}

	@Test
	public void testInfiniteRootKeepsFiniteChildrenAdjustable() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2, 3, 4));
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		final EdgeQueue<Integer> queue = EdgeQueue.create(1, partition);
		queue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), 3.0));
		queue.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), Double.POSITIVE_INFINITY));
		queue.addEdge(ExclusiveEdge.of(Edge.from(4).to(1), 1.0));
		final Edge<Integer> replace = Edge.from(2).to(1);
		final int component = partition.merge(1, 2);
		final EdgeQueue<Integer> merged = queues.merge(component, ImmutableList.of(
				Pair.of(queue, weighted(replace, 2.0))));

		final ExclusiveEdge<Integer> infinite = merged.popBestEdge().get();
		assertEquals(Double.POSITIVE_INFINITY, infinite.weight, 0.0);
		assertEquals(ImmutableList.of(replace), infinite.excluded);
		final ExclusiveEdge<Integer> next = merged.popBestEdge().get();
		assertEquals(1.0, next.weight, 0.0);
		assertEquals(ImmutableList.of(replace), next.excluded);
		assertEquals(-1.0, merged.popBestEdge().get().weight, 0.0);
	}

	@Test
	public void testNewEdgesDoNotInheritPendingContractionAdjustments() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2, 3, 4));
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		final EdgeQueue<Integer> original = EdgeQueue.create(1, partition);
		original.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), 30.0));
		original.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), 20.0));
		final Edge<Integer> replace = Edge.from(2).to(1);
		final int component = partition.merge(1, 2);
		final EdgeQueue<Integer> merged = queues.merge(component, ImmutableList.of(
				Pair.of(original, weighted(replace, 5.0))));
		merged.addEdge(ExclusiveEdge.of(Edge.from(4).to(2), 22.0));
		assertEquals(25.0, merged.popBestEdge().get().weight, 0.0);
		final ExclusiveEdge<Integer> fresh = merged.popBestEdge().get();
		assertEquals(22.0, fresh.weight, 0.0);
		assertTrue(fresh.excluded.isEmpty());
		final ExclusiveEdge<Integer> old = merged.popBestEdge().get();
		assertEquals(15.0, old.weight, 0.0);
		assertEquals(ImmutableList.of(replace), old.excluded);
	}

	@Test
	public void testLargeSequentialAdjustmentsStayFinite() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2, 3, 4));
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		EdgeQueue<Integer> queue = EdgeQueue.create(1, partition);
		queue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), -1.0e308));
		queue.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), -1.5e308));
		for (int node : new int[] {2, 4}) {
			final int component = partition.merge(1, node);
			queue = queues.merge(component, ImmutableList.of(
					Pair.of(queue, weighted(Edge.from(node).to(1), -1.0e308))));
		}
		assertEquals(1.0e308, queue.popBestEdge().get().weight, 0.0);
		assertEquals((-1.5e308 + 1.0e308) + 1.0e308, queue.popBestEdge().get().weight, 0.0);
		assertFalse(queue.popBestEdge().isPresent());
	}

	@Test
	public void testMixedScaleAdjustmentsPreserveHeapOrder() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2, 3, 4));
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		EdgeQueue<Integer> queue = EdgeQueue.create(1, partition);
		queue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), -1e16));
		queue.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), -1e16));
		// Adding 1 twice to -1e16 rounds back to -1e16, but adding 2 at once does not.
		for (int node : new int[] {2, 4}) {
			final int component = partition.merge(1, node);
			queue = queues.merge(component, ImmutableList.of(
					Pair.of(queue, weighted(Edge.from(node).to(1), -1.0))));
		}
		final double first = queue.popBestEdge().get().weight;
		assertTrue(first >= queue.popBestEdge().get().weight);
	}

	@Test
	public void testRandomMixedScaleMergesDrainInOrder() {
		final Random random = new Random(7324);
		for (int trial = 0; trial < 200; trial++) {
			final List<Integer> nodes = new ArrayList<>();
			for (int node = 0; node < 80; node++) nodes.add(node);
			final Partition<Integer> partition = Partition.singletons(nodes);
			final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
			EdgeQueue<Integer> queue = EdgeQueue.create(1, partition);
			for (int source = 40; source < 80; source++) {
				queue.addEdge(ExclusiveEdge.of(Edge.from(source).to(1), mixedScale(random)));
			}
			for (int node = 2; node < 40; node++) {
				final int component = partition.merge(1, node);
				queue = queues.merge(component, ImmutableList.of(
						Pair.of(queue, weighted(Edge.from(node).to(1), mixedScale(random)))));
				queue.addEdge(ExclusiveEdge.of(Edge.from(40 + node).to(node), mixedScale(random)));
			}
			double previous = Double.POSITIVE_INFINITY;
			int count = 0;
			for (var next = queue.popBestEdge(); next.isPresent(); next = queue.popBestEdge()) {
				assertTrue("Out of order in trial " + trial, next.get().weight <= previous);
				previous = next.get().weight;
				count++;
			}
			assertEquals(78, count);
		}
	}

	private static double mixedScale(Random random) {
		return Math.scalb(random.nextDouble() - 0.5, random.nextInt(600) - 300);
	}

	@Test(expected = ArithmeticException.class)
	public void testIncomparableWeightRangeIsRejected() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2));
		final EdgeQueue<Integer> queue = EdgeQueue.create(1, partition);
		queue.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), Double.MAX_VALUE));
		// Its gap below the other edge is not representable as a finite double.
		queue.addEdge(ExclusiveEdge.of(Edge.from(2).to(1), -Double.MAX_VALUE));
	}

	@Test
	public void testTiesFavorBestTreeAfterMergeAndRetainOtherCandidates() {
		final Partition<Integer> partition = Partition.singletons(ImmutableList.of(0, 1, 2, 3, 4));
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		final EdgeQueue<Integer> first = EdgeQueue.create(1, partition);
		final EdgeQueue<Integer> second = EdgeQueue.create(2, partition);
		first.addEdge(ExclusiveEdge.of(Edge.from(0).to(1), 8.0));
		first.addEdge(ExclusiveEdge.of(Edge.from(3).to(1), 8.0));
		second.addEdge(ExclusiveEdge.of(Edge.from(4).to(2), 10.0));
		final int component = partition.merge(1, 2);
		final EdgeQueue<Integer> merged = queues.merge(component, ImmutableList.of(
				Pair.of(first, weighted(Edge.from(2).to(1), 3.0)),
				Pair.of(second, weighted(Edge.from(1).to(2), 5.0))));
		final Arborescence<Integer> preferred = Arborescence.of(ImmutableMap.of(2, 4));
		final ExclusiveEdge<Integer> selected = merged.popBestEdge(preferred).get();
		assertEquals(Edge.from(4).to(2), selected.edge);
		assertEquals(5.0, selected.weight, 0.0);
		assertEquals(5.0, merged.popBestEdge().get().weight, 0.0);
		assertEquals(5.0, merged.popBestEdge().get().weight, 0.0);
		assertFalse(merged.popBestEdge().isPresent());
	}

	@Test
	public void testNestedMergesAgainstEagerWeightAndHistoryOracle() {
		final Random random = new Random(314159);
		final List<Integer> nodes = new ArrayList<>();
		for (int node = 0; node < 80; node++) nodes.add(node);
		final Partition<Integer> partition = Partition.singletons(nodes);
		final EdgeQueueMap<Integer> queues = new EdgeQueueMap<>(partition);
		final Map<Edge<Integer>, ExclusiveEdge<Integer>> expected = new HashMap<>();
		final List<Pair<EdgeQueue<Integer>, Integer>> pending = new ArrayList<>();
		for (int destination = 1; destination < 40; destination++) {
			final EdgeQueue<Integer> queue = EdgeQueue.create(destination, partition);
			for (int source = 40; source < nodes.size(); source++) {
				final ExclusiveEdge<Integer> edge = ExclusiveEdge.of(
						Edge.from(nodes.get(source)).to(nodes.get(destination)),
						(random.nextInt(2001) - 1000) / 8.0);
				queue.addEdge(edge);
				expected.put(edge.edge, edge);
			}
			pending.add(Pair.of(queue, destination));
		}
		while (pending.size() > 1) {
			final var first = pending.remove(random.nextInt(pending.size()));
			final var second = pending.remove(random.nextInt(pending.size()));
			final Edge<Integer> intoFirst = Edge.from(second.second).to(first.second);
			final Edge<Integer> intoSecond = Edge.from(first.second).to(second.second);
			final double firstWeight = (random.nextInt(201) - 100) / 8.0;
			final double secondWeight = (random.nextInt(201) - 100) / 8.0;
			for (Map.Entry<Edge<Integer>, ExclusiveEdge<Integer>> entry : expected.entrySet()) {
				final int destination = partition.componentOf(entry.getKey().destination);
				if (destination != first.second && destination != second.second) continue;
				final ExclusiveEdge<Integer> edge = entry.getValue();
				final List<Edge<Integer>> excluded = new ArrayList<>(edge.excluded);
				excluded.add(destination == first.second ? intoFirst : intoSecond);
				entry.setValue(ExclusiveEdge.of(edge.edge, excluded,
						edge.weight - (destination == first.second ? firstWeight : secondWeight)));
			}
			final int component = partition.merge(first.second, second.second);
			final EdgeQueue<Integer> merged = queues.merge(component, ImmutableList.of(
					Pair.of(first.first, weighted(intoFirst, firstWeight)),
					Pair.of(second.first, weighted(intoSecond, secondWeight))));
			pending.add(Pair.of(merged, component));
		}
		final EdgeQueue<Integer> result = pending.get(0).first;
		double previousWeight = Double.POSITIVE_INFINITY;
		while (!expected.isEmpty()) {
			final ExclusiveEdge<Integer> actual = result.popBestEdge().get();
			assertTrue(actual.weight <= previousWeight);
			previousWeight = actual.weight;
			final ExclusiveEdge<Integer> oracle = expected.remove(actual.edge);
			assertNotNull(oracle);
			assertEquals(oracle.weight, actual.weight, 0.0);
			assertEquals(oracle.excluded, actual.excluded);
		}
		assertFalse(result.popBestEdge().isPresent());
	}

	@Test
	public void testDeepExclusionHistoriesAreSharedAndTraversedWithoutRecursion() {
		final Edge<Integer> excluded = Edge.from(0).to(1);
		ExclusiveEdge.Exclusions<Integer> history = new ExclusiveEdge.Exclusions<>(ImmutableList.of());
		final ExclusiveEdge.Exclusions<Integer> singleton =
				new ExclusiveEdge.Exclusions<>(ImmutableList.of(excluded));
		for (int index = 0; index < 20000; index++) history = history.concat(singleton);
		assertEquals(20000, history.size());
		assertEquals(excluded, history.get(0));
		assertEquals(excluded, history.get(history.size() - 1));
		int count = 0;
		for (Edge<Integer> edge : history) {
			assertEquals(excluded, edge);
			count++;
		}
		assertEquals(20000, count);
	}
}
