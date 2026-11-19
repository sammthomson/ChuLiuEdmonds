package edu.cmu.cs.ark.cle.util;

/** Numeric policy used by the double specialization. */
public final class WeightMath {
	private WeightMath() {}

	public static boolean present(double weight) { return weight != Double.NEGATIVE_INFINITY; }
	public static double missing() { return Double.NEGATIVE_INFINITY; }
	public static double addScore(double a, double b) { return a + b; }

	public static double subtract(double a, double b) { return checked(a, b, a - b); }
	public static double combine(double a, double b) { return checked(a, b, a + b); }

	private static double checked(double a, double b, double result) {
		if (Double.isNaN(result) || (!Double.isFinite(result) && Double.isFinite(a) && Double.isFinite(b))) {
			throw new ArithmeticException("Unrepresentable contraction cost for " + a + " and " + b);
		}
		return result;
	}

	// Prediction and reconstruction can add rounded costs in different orders.
	public static final boolean EXACT = false;
}
