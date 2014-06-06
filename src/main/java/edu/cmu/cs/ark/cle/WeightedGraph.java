package edu.cmu.cs.ark.cle;

import java.util.Collection;

/**
 * @author sthomson@cs.cmu.edu
 */
public interface WeightedGraph<V> {
	public Collection<V> getNodes();

	public double getWeightOf(V source, V dest);
}
