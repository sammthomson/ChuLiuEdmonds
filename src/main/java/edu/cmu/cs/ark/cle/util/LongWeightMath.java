package edu.cmu.cs.ark.cle.util;

import java.util.NoSuchElementException;

/** Numeric policy used by the long specialization. No value is reserved for absent edges. */
public final class LongWeightMath {
	private LongWeightMath() {}

	public static boolean present(long weight) { return true; }
	public static long missing() { throw new NoSuchElementException("No such edge"); }
	public static long addScore(long a, long b) { return Math.addExact(a, b); }
	public static long subtract(long a, long b) { return Math.subtractExact(a, b); }
	public static long combine(long a, long b) { return Math.addExact(a, b); }
	public static final boolean EXACT = true;
}
