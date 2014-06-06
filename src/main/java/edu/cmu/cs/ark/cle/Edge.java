package edu.cmu.cs.ark.cle;

import com.google.common.base.Objects;

/** An edge in a directed graph. */
public class Edge<V> {
	public final V source;
	public final V destination;

	public Edge(V source, V destination) {
		this.source = source;
		this.destination = destination;
	}

	public static class EdgeBuilder<V> {
		public final V source;

		private EdgeBuilder(V source) {
			this.source = source;
		}

		public Edge<V> to(V destination) {
			return new Edge<V>(source, destination);
		}
	}

	public static <T> EdgeBuilder<T> from(T source) {
		return new EdgeBuilder<T>(source);
	}

	@Override public int hashCode() {
		return Objects.hashCode(source, destination);
	}

	@Override public String toString() {
		return Objects.toStringHelper(this)
				.add("source", source)
				.add("destination", destination).toString();
	}

	@Override public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Edge other = (Edge) obj;
		return this.source == other.source && this.destination == other.destination;
	}
}
