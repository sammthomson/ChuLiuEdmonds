package edu.cmu.cs.lti.ark.cle;

import com.google.common.primitives.Doubles;

import java.util.List;

/**
 * An edge, together with a list of edges that can't be in the final answer if 'edge' is.
 *
 * @author sthomson@cs.cmu.edu
 */
public class ExclusiveEdge implements Comparable<ExclusiveEdge> {
	public final Edge edge;
	public final List<Edge> excluded;
	public final double weight;

	public ExclusiveEdge(Edge edge, List<Edge> excluded, double weight) {
		this.edge = edge;
		this.excluded = excluded;
		this.weight = weight;
	}

	@Override public int compareTo(ExclusiveEdge exclusiveEdge) {
		return Doubles.compare(exclusiveEdge.weight, weight);
	}
}
