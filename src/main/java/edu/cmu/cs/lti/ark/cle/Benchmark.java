package edu.cmu.cs.lti.ark.cle;

import java.util.Random;

/**
 * @author sthomson@cs.cmu.edu
 */
public class Benchmark {
	public static void main(String[] args) {
		Random random = new Random(53242452);

		int burnin = 1000;
		int maxGraphSize = 10000;
		int numRunsPerSize = 100;

		for (int j = 0; j < burnin; j++) {
			final double[][] weights = randomGraph(200, random);
			ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
		}

		for (double n = 1; n < maxGraphSize; n *= 1.01) {
			final int n1 = new Double(n).intValue();
			final double[][] weights = randomGraph(n1, random);
			final long startTime = System.currentTimeMillis();
			for (int j = 0; j < numRunsPerSize; j++) {
				ChuLiuEdmonds.getMaxSpanningTree(weights, 0);
			}
			final long endTime = System.currentTimeMillis();
			System.out.println(n1 + "\t" + (endTime - startTime));
		}
	}

	public static double[][] randomGraph(int n, Random random) {
		double[][] result = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				result[i][j] = random.nextDouble();
			}
		}
		return result;
	}
}
