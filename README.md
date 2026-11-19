ChuLiuEdmonds
=============

Java implementation of Chu-Liu-Edmonds for maximum directed spanning arborescences,
including k-best enumeration. Negate edge weights to find a minimum arborescence.

## Build and test

Install [mise](https://mise.jdx.dev/), then run:

```sh
mise trust mise.toml
mise install java gradle
mise run test
mise run build
```

The repository pins Java 21 and Gradle. The build produces the library,
sources, and an `-all.jar` with runtime dependencies in `build/libs/`.
To run only the arborescence tests:

```sh
mise run test -- --tests edu.cmu.cs.ark.cle.ChuLiuEdmondsTest --tests edu.cmu.cs.ark.cle.KBestArborescencesTest
```

The regression tests independently enumerate small rooted trees to check
optimum scores, required/banned edges, and complete k-best rankings.

Run `mise run benchmark` for a deliberately contraction-heavy sparse graph family.
It reports the median of three runs after a small warmup; graph construction is
outside the timed region. This is a focused diagnostic, not a general benchmark.

## Choosing a weight type

The existing `ChuLiuEdmonds` and `KBestArborescences` APIs use `double`.
For exact integer comparisons and scores, use their `Long`-prefixed counterparts:

```java
var graph = LongSparseWeightedGraph.from(List.of(
    LongWeighted.weighted(Edge.from(0).to(1), 9_007_199_254_740_993L)));
long score = LongChuLiuEdmonds.getMaxArborescence(graph, 0).weight;
var ranked = LongKBestArborescences.getKBestArborescences(graph, 0, 3);
```

Both use primitive weights throughout the heap. Gradle generates the `long` classes
from the same algorithm sources; neither API converts inputs to the other type.
The `long` implementation throws `ArithmeticException` if an intermediate cost,
score, or ranking difference overflows. Every `long`, including `Long.MIN_VALUE`,
is a valid weight. Use a sparse graph to omit edges; every dense matrix entry
represents an edge. Looking up an absent long edge throws `NoSuchElementException`.

With `double`, ordinary queue insertions and removals preserve input priorities
when no contraction costs apply. During contractions, combined costs may round
differently from successive subtractions. On mixed-scale inputs this can change
the selected tree, the reported score, and the k-best ordering. An undefined or
overflowing contraction cost throws `ArithmeticException`.

## Algorithm

The implementation combines Tarjan's cycle contractions with the **growth path**
processing order of Gabow, Galil, Spencer, and Tarjan. After selecting an incoming
edge, it continues at that edge's source. Revisiting a component on the current
path identifies a cycle, whose path suffix is contracted. A path ends at a
component with no incoming edge or at an already completed path.

This removes the separate weak-component partition and incoming-edge map.
The contraction partition and exclusion-based reconstruction are retained;
this is the **Tarjan+Path** variant, not the full GGST algorithm with its
additional data structures and improved asymptotic bound.

Incoming-edge queues use a leftist heap. Contraction subtracts the replaced edge's
weight from the queue root and defers the same cost for its children. The heap pushes
pending costs before changing children, so contraction can meld the old queues
instead of scanning and reinserting their edges.
Exclusion histories are immutable concatenation trees shared between edges;
they are traversed iteratively when reconstructing the answer. Edges that become
internal to a contracted component are discarded lazily at the queue head.

This removes bulk per-edge contraction updates without adding GGST's exit lists,
passive lists, or specialized active forest. It does not claim the full GGST
bound: tie handling and exclusion-based reconstruction can still do extra work.

### Why defer contraction costs?

We use lazy offsets because nested contractions repeatedly adjust the same queued
edges. Eager subtraction visits all of them at each contraction, including edges
that will never be selected. To isolate that cost, we compared eager and lazy
subtraction with the same leftist heap, traversal, and exclusion handling.

| Graph (nodes) | `double` lazy (ms) | `double` eager (ms) | `long` lazy (ms) | `long` eager (ms) |
|---|---:|---:|---:|---:|
| Sparse DAG (5,000) | 15.25 | 13.01 | 14.59 | 15.75 |
| Sparse cyclic (5,000) | 14.03 | 12.18 | 12.14 | 12.78 |
| Dense cyclic (250) | 13.26 | 13.47 | 11.10 | 11.94 |
| Nested contractions (1,001) | 2.57 | 3.44 | 1.84 | 3.56 |
| Nested contractions (2,001) | 5.17 | 10.96 | 3.80 | 11.34 |

These are median times for a rooted maximum-tree query, excluding graph
construction: three fresh Java 21 JVMs per variant, five warmups and seven timed
runs per workload, then the median of the three run medians. Lazy offsets were
2.1 times faster for `double` and 3.0 times faster for `long` on the larger nested
graph. Neither strategy won consistently on the other graph families.

Eager updates also showed no advantage in **tree selection** on the sampled
inputs. We enumerated every rooted tree on 2,000 five-node graphs per weight
distribution and summed the original weights exactly. On mixed-scale integers
(small integers mixed with multiples of 2^52), each `double` strategy selected
a suboptimal tree in 13 of 2,000 cases; both `long` strategies returned an
optimum with its exact score in every case. On random fractions in `[-1, 1)`,
both `double` strategies always selected an optimum, but each reported a score
different from the correctly rounded sum of its selected edges in 284 of 2,000
cases. These counts describe the samples, not a general error bound: eager
subtraction changes the order of rounding, but does not make `double` exact.

### References

- R. E. Tarjan, "Finding Optimum Branchings," *Networks* 7 (1977), 25-35.
  [DOI](https://doi.org/10.1002/net.3230070103).
- H. N. Gabow, Z. Galil, T. Spencer, and R. E. Tarjan, "Efficient algorithms for
  finding minimum spanning trees in undirected and directed graphs,"
  *Combinatorica* 6 (1986), 109-122.
  [DOI](https://doi.org/10.1007/BF02579168).
- O. Kissig and C. Weyand, "Efficiently Computing Directed Minimum Spanning Trees,"
  Sections 2.3 and 3 explain the growth path and the Tarjan+Path variant.
  [arXiv:2208.02590](https://arxiv.org/abs/2208.02590).
