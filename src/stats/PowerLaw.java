package stats;

import java.util.Random;

/**
 * @deprecated essentially the same as Zipf, so use that instead
 */
public class PowerLaw
        // implements IndelModel
    {

    private double alpha; // Power law exponent
    private Random rand = null;
    private long seed = System.currentTimeMillis();
    private int maxK; // Maximum rank value

    /**
     * Define a Power Law distribution
     *
     * @param alpha the exponent parameter for the distribution
     * @param seed the random seed for sampling
     * @param maxK the maximum rank value for sampling
     */
    public PowerLaw(double alpha, long seed, int maxK) {
        this.alpha = alpha;
        this.seed = seed;
        this.maxK = maxK;
    }

    /**
     * Probability mass function (PMF)
     *
     * @param k the rank
     * @return the probability of rank k
     */
    public double p(int k) {
        return Math.pow(k, -alpha) / normalizationConstant();
    }

    /**
     * Cumulative distribution function (CDF)
     *
     * @param k the rank
     * @return the cumulative probability up to rank k
     */
    public double cdf(int k) {
        double sum = 0;
        for (int i = 1; i <= k; i++) {
            sum += p(i);
        }
        return sum;
    }

    /**
     * Sample from the Power Law distribution
     *
     * @return the sampled rank
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
     * Calculate the normalization constant
     *
     * @return the normalization constant
     */
    private double normalizationConstant() {
        double sum = 0;
        for (int i = 1; i <= maxK; i++) {
            sum += Math.pow(i, -alpha);
        }
        return sum;
    }

    public String getTrAVIS() {
        return "PowerLaw " + alpha;
    }

    public void setSeed(long seed) {
        this.seed = seed;
        this.rand = new Random(seed);
    }


    public static void main(String[] args) {
        // Example usage with different alpha and maxK values
        for (double alpha = 1.0; alpha < 5.5; alpha += 0.5) {
            PowerLaw powerLaw = new PowerLaw(alpha, 0, 1000); // Set maxK to 5000
            System.out.println("alpha = " + alpha);
            int[] cnt = new int[powerLaw.maxK + 1];
            for (int i = 0; i < 1000; i++) {
                cnt[powerLaw.sample()] += 1;
            }
            for (int j = 1; j <= 20; j++) {
                System.out.print(cnt[j] + "\t");
            }
            System.out.println();
        }
    }
}

