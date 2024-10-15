package stats;

import java.util.Random;

public class Zipf implements IndelModel {

    private double s; // Exponent parameter for the Zipf distribution
    private Random rand = null;
    private long seed = System.currentTimeMillis();
    private int maxK; // Maximum value for the sampling range

    /**
     * Defines a Zipf distribution
     *
     * @param s The exponent parameter of the distribution
     * @param seed The random seed for sampling
     * @param maxK The maximum value for the sampling range
     */
    public Zipf(double s, long seed, int maxK) {
        this.s = s;
        this.seed = seed;
        this.maxK = maxK;
    }

    public String toString() {
        return "Zipf(lambda=" + s + ", seed=" + seed + ", maxK=" + maxK + ")";
    }

    /**
     * Probability mass function (PMF)
     *
     * @param k The rank
     * @return The probability of rank k
     */
    public double p(int k) {
        return 1.0 / Math.pow(k, s) / zeta(s);
    }

    /**
     * Cumulative distribution function (CDF)
     *
     * @param k The rank up to which the cumulative probability is calculated
     * @return The cumulative probability of ranks up to and including k
     */
    public double cdf(int k) {
        double sum = 0;
        for (int i = 1; i <= k; i++) {
            sum += p(i);
        }
        return sum;
    }

    /**
     * Samples from the Zipf distribution
     *
     * @return A rank sampled from the distribution
     */
    public int sample() {
        if (rand == null) {
            rand = new Random(seed);
        }
        double toss = rand.nextDouble();
        double sum = 0;
        for (int i = 1; i <= maxK; i++) {
            sum += p(i);
            if (sum > toss) {
                return i;
            }
        }
        return maxK;
    }

    /**
     * Computes an approximation of the Riemann zeta function value
     *
     * @param s The exponent parameter
     * @return The approximate value of the zeta function
     */
    private double zeta(double s) {
        double sum = 0;
        for (int i = 1; i <= maxK; i++) {
            sum += 1.0 / Math.pow(i, s);
        }
        return sum;
    }

    public static void main(String[] args) {
        // Example: Using different s and maxK values
        for (double s = 1.0; s < 5.5; s += 0.5) {
            Zipf zipf = new Zipf(s, 1, 1000); // Set maxK to 5000
            System.out.println("s = " + s);
            int[] cnt = new int[zipf.maxK + 1];
            for (int i = 0; i < 1000; i++) {
                cnt[zipf.sample()] += 1;
            }
            for (int j = 1; j <= 20; j++) {
                System.out.print(cnt[j] + "\t");
            }
            System.out.println();
        }
    }
}
