package stats;

import java.util.Random;
import java.util.ArrayList;
import java.util.List;

public class PowerLawCon {

    private double alpha; // Power law exponent
    private Random rand = null;
    private long seed = System.currentTimeMillis();
    private double min; // Minimum value for the distribution
    private double max; // Maximum value for the distribution
    private double normalizationConstant;

    /**
     * Define a continuous Power Law distribution
     *
     * @param alpha the exponent parameter for the distribution
     * @param seed  the random seed for sampling
     * @param min   the minimum value for the distribution
     * @param max   the maximum value for the distribution
     */
    public PowerLawCon(double alpha, long seed, double min, double max) {
        if (min <= 0 || max <= 0 || min >= max) {
            throw new IllegalArgumentException("min and max must be positive, and min < max.");
        }
        this.alpha = alpha;
        this.seed = seed;
        this.min = min;
        this.max = max;
        this.normalizationConstant = calculateNormalizationConstant();
    }

    /**
     * Probability density function (PDF)
     *
     * @param x the value
     * @return the probability density at x
     */
    public double pdf(double x) {
        if (x < min || x > max) {
            return 0.0;
        }
        return (alpha - 1) * Math.pow(x, -alpha) / normalizationConstant;
    }

    /**
     * Cumulative distribution function (CDF)
     *
     * @param x the value
     * @return the cumulative probability up to x
     */
    public double cdf(double x) {
        if (x < min) {
            return 0.0;
        }
        if (x > max) {
            return 1.0;
        }
        double cumulative = (Math.pow(x, 1 - alpha) - Math.pow(min, 1 - alpha)) /
                (Math.pow(max, 1 - alpha) - Math.pow(min, 1 - alpha));
        return cumulative;
    }

    /**
     * Sample from the continuous Power Law distribution
     *
     * @return the sampled value
     */
    public double sample() {
        if (rand == null) {
            rand = new Random(seed);
        }
        double u = rand.nextDouble();
        double exponent = 1 - alpha;
        return Math.pow((u * (Math.pow(max, exponent) - Math.pow(min, exponent)) + Math.pow(min, exponent)), 1 / exponent);
    }

    /**
     * Calculate the normalization constant
     *
     * @return the normalization constant
     */
    private double calculateNormalizationConstant() {
        return Math.pow(max, 1 - alpha) - Math.pow(min, 1 - alpha);
    }

    public String getTrAVIS() {
        return "PowerLawCon " + alpha + " [" + min + ", " + max + "]";
    }

    /**
     * Estimate alpha (exponent) using Maximum Likelihood Estimation (MLE)
     *
     * @param samples a list of sampled data
     * @param min     the minimum value of the distribution
     * @return the estimated alpha
     */
    public static double estimateAlpha(List<Double> samples, double min) {
        if (samples == null || samples.isEmpty()) {
            throw new IllegalArgumentException("Samples list cannot be null or empty.");
        }

        // Calculate the log of each sample
        double logSum = 0.0;
        for (double sample : samples) {
            if (sample < min) {
                throw new IllegalArgumentException("All samples must be >= min.");
            }
            logSum += Math.log(sample / min);
        }

        // Estimate alpha
        double n = samples.size();
        return 1 + n / logSum;
    }

    /**
     * Compute min and max dynamically from the samples
     *
     * @param samples the list of sampled data
     * @return an array where [0] is min and [1] is max
     */
    public static double[] computeMinMax(List<Double> samples) {
        if (samples == null || samples.isEmpty()) {
            throw new IllegalArgumentException("Samples list cannot be null or empty.");
        }

        double min = samples.stream().min(Double::compare).orElse(0.0);
        double max = samples.stream().max(Double::compare).orElse(0.0);

        return new double[]{min, max};
    }

    public static void main(String[] args) {
        // Generate samples
        double alpha = 1.5;
        List<Double> samples = new ArrayList<>();
        PowerLawCon powerLawCon = new PowerLawCon(alpha, 0, 0.01, 100.0);

        for (int i = 0; i < 1000; i++) {
            samples.add(powerLawCon.sample());
        }

        // Compute min and max using the method
        double[] minMax = computeMinMax(samples);
        double min = minMax[0];
        double max = minMax[1];

        System.out.println("Min: " + min);
        System.out.println("Max: " + max);

        // Estimate alpha using the computed min
        double estimatedAlpha = estimateAlpha(samples, min);
        System.out.println("Estimated alpha: " + estimatedAlpha);
    }
}

