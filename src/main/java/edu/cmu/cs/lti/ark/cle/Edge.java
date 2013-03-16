package edu.cmu.cs.lti.ark.cle;

import com.google.common.base.Objects;

/** An edge in a directed graph. */
public class Edge {
	public final int source;
	public final int destination;

	public Edge(int source, int destination) {
		this.source = source;
		this.destination = destination;
	}

	@Override public int hashCode() {
		return Objects.hashCode(source, destination);
	}

	@Override public String toString() {
		return Objects.toStringHelper(this)
				.add("source", source)
				.add("destination", destination).toString();
	}
}