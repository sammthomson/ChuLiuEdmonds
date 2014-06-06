package edu.cmu.cs.ark.cle;

import com.google.common.collect.ImmutableMap;

import java.util.Map;

/**
 * @author sthomson@cs.cmu.edu
 */
public class Arborescence {
	final static Arborescence EMPTY = new Arborescence(ImmutableMap.<Integer, Integer>of());

	final Map<Integer,Integer> parents;

	public Arborescence(Map<Integer, Integer> parents) {
		this.parents = parents;
	}
}
