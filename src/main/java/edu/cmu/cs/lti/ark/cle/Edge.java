package edu.cmu.cs.lti.ark.cle;

import com.google.common.base.Objects;

/** An edge in a directed graph. */
public class Edge<V> {
	public final V source;
	public final V destination;

	public Edge(V source, V destination) {
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