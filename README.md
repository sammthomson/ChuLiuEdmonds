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

Incoming-edge queues use a small leftist heap. Each root stores its absolute weight,
and each other node stores its nonnegative gap below its parent, so adjusting a root
adjusts its whole subtree. Floating-point rounding can then never place a child above
its parent. Finite children of a `+Infinity` parent keep absolute weights instead.
Finite weights whose gaps overflow a `double`, and adjustments yielding `∞ - ∞`,
are rejected with an `ArithmeticException`.
Contraction melds the old queues instead of scanning and reinserting their edges.
Exclusion histories are immutable concatenation trees shared between edges;
they are traversed iteratively when reconstructing the answer. Edges that become
internal to a contracted component are discarded lazily at the queue head.

This removes bulk per-edge contraction updates without adding GGST's exit lists,
passive lists, or specialized active forest. It does not claim the full GGST
bound: tie handling and exclusion-based reconstruction can still do extra work.

References:

- R. E. Tarjan, "Finding Optimum Branchings," *Networks* 7 (1977), 25-35.
  [DOI](https://doi.org/10.1002/net.3230070103).
- H. N. Gabow, Z. Galil, T. Spencer, and R. E. Tarjan, "Efficient algorithms for
  finding minimum spanning trees in undirected and directed graphs,"
  *Combinatorica* 6 (1986), 109-122.
  [DOI](https://doi.org/10.1007/BF02579168).
- O. Kissig and C. Weyand, "Efficiently Computing Directed Minimum Spanning Trees,"
  Sections 2.3 and 3 explain the growth path and the Tarjan+Path variant.
  [arXiv:2208.02590](https://arxiv.org/abs/2208.02590).
